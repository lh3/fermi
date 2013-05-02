#include <assert.h>
#include <string.h>
#include <pthread.h>
#include "rld.h"
#include "fermi.h"
#include "kvec.h"

#define SUF_LEN 4

static fmintv_t descend(const rld_t *e, int suf_len, int suf)
{
	int i;
	fmintv_t ok[6], ik;
	fm6_set_intv(e, (suf&3) + 1, ik);
	for (i = 1; i < suf_len; ++i) {
		fm6_extend(e, &ik, ok, 1);
		ik = ok[(suf>>i*2&3) + 1];
	}
	return ik;
}

static void collect_tips(const rld_t *e, uint64_t *sub, const fmintv_t *_ik, fmintv_v *stack)
{
	stack->n = 0;
	kv_push(fmintv_t, *stack, *_ik);
	while (stack->n) {
		uint64_t k, *p, x;
		fmintv_t ik, ok[6];
		int c;
		ik = kv_pop(*stack);
		fm6_extend(e, &ik, ok, 1);
		if (ok[0].x[2]) {
			for (k = 0; k < ok[0].x[2]; ++k) {
				x = e->mcnt[1] - 1 - (k + ok[0].x[0]);
				p = sub + (x>>6);
				x = 1ULL<<(x&0x3f);
				__sync_or_and_fetch(p, x);
			}
		}
		for (c = 1; c <= 4; ++c)
			if (ok[c].x[2]) kv_push(fmintv_t, *stack, ok[c]);
	}
}

static void contrast_core(const rld_t *eref, const rld_t *eqry, uint64_t *sub, int kmer, int min_occ, int suf_len, int suf)
{
	fmintv_v stack[2], tstack;
	fmintv_t ik[2], ok[2][6];
	const rld_t *e[2];
	int i, c;

	e[0] = eref; e[1] = eqry;
	kv_init(tstack); // temporary stack
	for (i = 0; i < 2; ++i) {
		kv_init(stack[i]);
		ik[i] = descend(e[i], suf_len, suf);
		ik[i].info = suf_len;
	}
	if (ik[0].x[2] == 0 || ik[1].x[2] < min_occ) return;
	kv_push(fmintv_t, stack[0], ik[0]);
	kv_push(fmintv_t, stack[1], ik[1]);
	while (stack[0].n) { // stack[0] and stack[1] are always of the same size
		ik[0] = kv_pop(stack[0]);
		ik[1] = kv_pop(stack[1]); // it is always true that ik[1].x[2] >= min_occ
		if (ik[0].x[2] == 0) collect_tips(e[1], sub, &ik[1], &tstack);
		else if (ik[1].info >= kmer) continue;
		else {
			fm6_extend(e[0], &ik[0], ok[0], 1);
			fm6_extend(e[1], &ik[1], ok[1], 1);
			for (c = 1; c <= 4; ++c) {
				if (ok[1][c].x[2] < min_occ) continue;
				ok[1][c].info = ik[1].info + 1;
				kv_push(fmintv_t, stack[0], ok[0][c]);
				kv_push(fmintv_t, stack[1], ok[1][c]);
			}
		}
	}
	free(stack[0].a); free(stack[1].a); free(tstack.a);
}

typedef struct {
	int tid, n_suf, k, min_occ;
	uint32_t *suf;
	const rld_t *eref, *eqry;
	uint64_t *sub;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int i;
	for (i = 0; i < w->n_suf; ++i)
		contrast_core(w->eref, w->eqry, w->sub, w->k, w->min_occ, SUF_LEN, w->suf[i]);
	return 0;
}

uint64_t *fm6_contrast(const rld_t *eref, const rld_t *eqry, int k, int min_occ, int n_threads)
{
	worker_t *w;
	pthread_t *tid;
	pthread_attr_t attr;
	int i, j, max_suf;
	uint64_t *sub;

	assert(k > SUF_LEN);
	if (!(n_threads&1)) --n_threads; // make it to an odd number
	if (n_threads == 0) n_threads = 1;
	tid = malloc(n_threads * sizeof(pthread_t));
	sub = calloc((eqry->mcnt[1] + 63) / 64, 8);
	w = calloc(n_threads, sizeof(worker_t));
	max_suf = ((1<<SUF_LEN*2) + n_threads - 1) / n_threads;
	for (j = 0; j < n_threads; ++j) {
		w[j].suf = calloc(max_suf, 4);
		w[j].k = k, w[j].min_occ = min_occ, w[j].tid = j;
		w[j].eref = eref, w[j].eqry = eqry, w[j].sub = sub;
	}
	max_suf = 1<<SUF_LEN*2;
	for (i = 0, j = 0; i < max_suf; ++i) { // assign seqs and hash tables
		w[j].suf[w[j].n_suf++] = i;
		if (++j == n_threads) j = 0;
	}
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	for (j = 0; j < n_threads; ++j) free(w[j].suf);
	free(w); free(tid);
	return sub;
}
