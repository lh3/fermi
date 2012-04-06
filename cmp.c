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
				x = k + ok[0].x[0];
				p = sub + (x>>6);
				x = 1ULL<<(x&0x3f);
				__sync_or_and_fetch(p, x);
			}
		}
		for (c = 1; c <= 4; ++c)
			if (ok[c].x[2]) kv_push(fmintv_t, *stack, ok[c]);
	}
}

static void contrast_core(const rld_t *e[2], uint64_t *sub[2], int kmer, int min_occ, int suf_len, int suf)
{
	fmintv_v stack[2], tstack;
	fmintv_t ik[2], ok[2][6];
	int i, c;

	kv_init(tstack);
	for (i = 0; i < 2; ++i) {
		kv_init(stack[i]);
		ik[i] = descend(e[i], suf_len, suf);
		ik[i].info = suf_len;
		kv_push(fmintv_t, stack[i], ik[i]);
	}
	while (stack[0].n) {
		ik[0] = kv_pop(stack[0]);
		ik[1] = kv_pop(stack[1]);
		if (ik[0].x[2] == 0)      collect_tips(e[1], sub[1], &ik[1], &tstack);
		else if (ik[1].x[2] == 0) collect_tips(e[0], sub[0], &ik[0], &tstack);
		else if (ik[0].info >= kmer) continue;
		else {
			fm6_extend(e[0], &ik[0], ok[0], 1);
			fm6_extend(e[1], &ik[1], ok[1], 1);
			for (c = 1; c <= 4; ++c) {
				if (ok[0][c].x[2] < min_occ && ok[1][c].x[2] < min_occ) continue;
				ok[0][c].info = ik[0].info + 1;
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
	const rld_t *e[2];
	uint64_t *sub[2];
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int i;
	for (i = 0; i < w->n_suf; ++i)
		contrast_core(w->e, w->sub, w->k, w->min_occ, SUF_LEN, w->suf[i]);
	return 0;
}

void fm6_contrast(rld_t *const e[2], int k, int min_occ, int n_threads, uint64_t *sub[2])
{
	worker_t *w;
	pthread_t *tid;
	pthread_attr_t attr;
	int i, j, max_suf;

	assert(k > SUF_LEN);
	if (!(n_threads&1)) --n_threads; // make it to an odd number
	if (n_threads == 0) n_threads = 1;
	tid = malloc(n_threads * sizeof(pthread_t));
	sub[0] = calloc((e[0]->mcnt[1] + 63) / 64, 8);
	sub[1] = calloc((e[1]->mcnt[1] + 63) / 64, 8);
	w = calloc(n_threads, sizeof(worker_t));
	max_suf = ((1<<SUF_LEN*2) + n_threads - 1) / n_threads;
	for (j = 0; j < n_threads; ++j) {
		w[j].suf = calloc(max_suf, 4);
		w[j].k = k, w[j].min_occ = min_occ, w[j].tid = j;
		w[j].e[0] = e[0], w[j].e[1] = e[1];
		w[j].sub[0] = sub[0], w[j].sub[1] = sub[1];
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
}

int64_t fm6_sub_conv(int64_t n_seqs, uint64_t *sub, const uint64_t *rank)
{
	uint64_t i, k, *tmp, n_sel;
	tmp = calloc((n_seqs + 63) / 64, 8);
	for (i = n_sel = 0; i < n_seqs; ++i) {
		if (sub[i>>6]>>(i&0x3f)&1) {
			k = rank[i]>>2;
			tmp[k>>6] |= 1ULL<<(k&0x3f);
			++n_sel;
		}
	}
	memcpy(sub, tmp, (n_seqs + 63) / 64 * 8);
	free(tmp);
	for (i = 0; i < n_seqs; i += 2)
		assert(((sub[i>>6]>>(i&0x3f) ^ sub[i>>6]>>((i^1)&0x3f)) & 1) == 0);
	return n_sel;
}
