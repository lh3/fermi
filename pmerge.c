#include <pthread.h>
#include <stdio.h>
#include "rld.h"
#include "utils.h"

#define BLOCK_SIZE 0x1000000

typedef struct {
	int step, n;
	int64_t size, i, k, x, *buf;
	const rld_t *e0, *e1;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	uint64_t *ok;
	if (w->x >= w->e1->mcnt[1]) return 0;
	ok = alloca(8 * w->e0->asize);
	w->n = 0;
	w->buf[w->n++] = w->i + w->k + 1;
	for (;;) {
		int c = rld_rank1a(w->e1, w->k, ok);
		if (c == 0) {
			w->x += w->step;
			if (w->x >= w->e1->mcnt[1]) break;
			w->k = w->x;
			w->i = w->e0->mcnt[1] - 1;
		} else {
			w->k = w->e1->cnt[c] + ok[c] - 1;
			rld_rank1a(w->e0, w->i, ok);
			w->i = w->e0->cnt[c] + ok[c] - 1;
		}
		if (w->n == w->size) return 0;
		else w->buf[w->n++] = w->k + w->i + 1;
	}
	return 0;
}

uint64_t *fm_compute_gap_bits(const rld_t *e0, const rld_t *e1, int n_threads)
{
	uint64_t *bits;
	int64_t rest = e1->mcnt[0];
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *data;
	int j;

	bits = xcalloc((e0->mcnt[0] + e1->mcnt[0] + 63) / 64, 8);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	data = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (j = 0; j < n_threads; ++j) {
		worker_t *w = data + j;
		w->e0 = e0; w->e1 = e1;
		w->step = n_threads;
		w->size = ((BLOCK_SIZE < rest? BLOCK_SIZE : rest) + n_threads - 1) / n_threads;
		w->k = w->x = j;
		w->i = w->e0->mcnt[1] - 1;
		w->buf = malloc(w->size * 8);
	}
	while (rest) {
		for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, data + j);
		for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
		for (j = 0; j < n_threads; ++j) {
			int i;
			worker_t *w = data + j;
			for (i = 0; i < w->n; ++i)
				bits[w->buf[i]>>6] |= 1ull<<(w->buf[i]&0x3f);
			rest -= w->n;
		}
	}
	free(data); free(tid);
	return bits;
}
