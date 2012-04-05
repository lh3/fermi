#include <assert.h>
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

static void contrast_core(const rld_t *ref, const rld_t *src, uint64_t *set, int kmer, int min_occ, int suf_len, int suf)
{
	fmintv_v rstack, sstack, tstack;
	fmintv_t rik, sik, rok[6], sok[6];

	kv_init(rstack); kv_init(sstack); kv_init(tstack);
	rik = descend(ref, suf_len, suf);
	sik = descend(src, suf_len, suf); sik.info = suf_len;
	kv_push(fmintv_t, rstack, rik);
	kv_push(fmintv_t, sstack, sik);
	while (sstack.n) {
		int c;
		sik = kv_pop(sstack);
		rik = kv_pop(rstack); // rstack and sstack always has the same number of elements
		if (rik.x[2] == 0) {
			tstack.n = 0;
			kv_push(fmintv_t, tstack, sik);
			while (tstack.n) {
				uint64_t k, *p, x;
				sik = kv_pop(tstack);
				fm6_extend(src, &sik, sok, 1);
				if (sok[0].x[2]) {
					for (k = 0; k < sok[0].x[2]; ++k) {
						x = k + sok[0].x[0];
						p = set + (x>>6);
						x = 1ULL<<(x&0x3f);
						__sync_or_and_fetch(p, x);
					}
				}
				for (c = 1; c <= 4; ++c)
					if (sok[c].x[2]) kv_push(fmintv_t, tstack, sok[c]);
			}
		} else if (sik.info >= kmer) {
			continue;
		} else {
			fm6_extend(src, &sik, sok, 1);
			fm6_extend(ref, &rik, rok, 1);
			for (c = 1; c <= 4; ++c) {
				if (sok[c].x[2] < min_occ) continue;
				sok[c].info = sik.info + 1;
				kv_push(fmintv_t, sstack, sok[c]);
				kv_push(fmintv_t, rstack, rok[c]);
			}
		}
	}
	free(rstack.a); free(sstack.a); free(tstack.a);
}

typedef struct {
	int tid, n_suf, k, min_occ;
	uint32_t *suf;
	const rld_t *ref, *src;
	uint64_t *set;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int i;
	for (i = 0; i < w->n_suf; ++i)
		contrast_core(w->ref, w->src, w->set, w->k, w->min_occ, SUF_LEN, w->suf[i]);
	return 0;
}

uint64_t *fm6_contrast(const rld_t *ref, const rld_t *src, int k, int min_occ, int n_threads)
{
	worker_t *w;
	pthread_t *tid;
	pthread_attr_t attr;
	uint64_t *set;
	int i, j, max_suf;

	assert(k > SUF_LEN);
	if (!(n_threads&1)) --n_threads; // make it to an odd number
	if (n_threads == 0) n_threads = 1;
	tid = malloc(n_threads * sizeof(pthread_t));
	set = calloc((src->mcnt[1] + 63) / 64, 8);
	w = calloc(n_threads, sizeof(worker_t));
	max_suf = ((1<<SUF_LEN*2) + n_threads - 1) / n_threads;
	for (j = 0; j < n_threads; ++j) {
		w[j].suf = calloc(max_suf, 4);
		w[j].ref = ref, w[j].src = src, w[j].tid = j;
		w[j].k = k, w[j].min_occ = min_occ;
		w[j].set = set;
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
	return set;
}
