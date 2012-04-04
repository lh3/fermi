#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include "priv.h"

typedef struct {
	const rld_t *e;
	uint64_t *sorted;
	int start, step;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int64_t i;
	kstring_t s;
	s.l = s.m = 0; s.s = 0;
	for (i = w->start<<1; i < w->e->mcnt[1]; i += w->step<<1) {
		fmintv_t k2;
		uint64_t k, l;
		int contained, flag;
		s.l = 0;
		k = fm6_retrieve(w->e, i, &s, &k2, &contained);
		flag = (contained != 0)<<1 | (k2.x[2] > 1 && k != k2.x[0]);
		w->sorted[k] = i<<2 | flag;
		assert(k >= k2.x[0] && k < k2.x[0] + k2.x[2]);
		if (k2.x[0] != k2.x[1]) { // seq and reverse complement are different
			for (l = 0; l < k2.x[2]; ++l)
				if (k == k2.x[0] + l) break;
			w->sorted[k2.x[1] + l] = (i|1)<<2 | flag;
		} else w->sorted[k+1] = (i|1)<<2 | flag;
	}
	free(s.s);
	return 0;
}

uint64_t *fm6_seqsort(const rld_t *e, int n_threads)
{
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	uint64_t *sorted;
	int j;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	sorted = calloc(e->mcnt[1], 8);
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e = e;
		ww->step = n_threads;
		ww->start = j;
		ww->sorted = sorted;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(w); free(tid);
	{
		uint64_t i, cnt0 = 0, n_contained = 0, n_dups = 0;
		for (i = 0; i < e->mcnt[1]; ++i)
			if (sorted[i] == 0) ++cnt0;
			else if (sorted[i]&2) ++n_contained;
			else if (sorted[i]&1) ++n_dups; // if contained, not counted as dups
		fprintf(stderr, "[M::%s] #zeros=%ld, #contained=%ld, #duplicates=%ld\n", __func__,
				(long)cnt0, (long)n_contained, (long)n_dups);
	}
	return sorted;
}
