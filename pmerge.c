#include <pthread.h>
#include <stdio.h>
#include "fermi.h"
#include "rld.h"
#include "utils.h"

#define BLOCK_SIZE 0x1000000
#define PWRITE_MIN 4

typedef struct {
	int step, n;
	int64_t size, i, k, x, *buf;
	const rld_t *e0, *e1;
} worker_t;

typedef struct {
	int tid, n_threads;
	worker_t *w;
	uint64_t *bits;
} writer_t;

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

static void *writer(void *data)
{
	writer_t *w = (writer_t*)data;
	int i, j;
	for (j = 0; j < w->n_threads; ++j) {
		worker_t *wo = w->w + j;
		for (i = 0; i < wo->n; ++i) {
			uint64_t x = wo->buf[i];
			if ((x>>6) % w->n_threads == w->tid)
				w->bits[x>>6] |= 1ull<<(x&0x3f);
		}
	}
	return 0;
}

uint64_t *fm_compute_gap_bits(const rld_t *e0, const rld_t *e1, int n_threads)
{
	uint64_t *bits;
	int64_t rest = e1->mcnt[0];
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *wo;
	writer_t *wr;
	int j;
	double tcpu, treal;

	bits = xcalloc((e0->mcnt[0] + e1->mcnt[0] + 63) / 64, 8);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	wo = (worker_t*)calloc(n_threads, sizeof(worker_t));
	wr = (writer_t*)calloc(n_threads, sizeof(writer_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (j = 0; j < n_threads; ++j) {
		worker_t *w = wo + j;
		writer_t *r = wr + j;
		w->e0 = e0; w->e1 = e1;
		w->step = n_threads;
		w->size = ((BLOCK_SIZE < rest? BLOCK_SIZE : rest) + n_threads - 1) / n_threads;
		w->k = w->x = j;
		w->i = w->e0->mcnt[1] - 1;
		w->buf = malloc(w->size * 8);
		r->n_threads = n_threads;
		r->tid = j;
		r->bits = bits;
		r->w = wo;
	}
	tcpu = cputime(); treal = realtime();
	while (rest) {
		for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, wo + j);
		for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
		if (n_threads >= PWRITE_MIN) {
			for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, writer, wr + j);
			for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
			for (j = 0; j < n_threads; ++j) rest -= wo[j].n;
		} else {
			for (j = 0; j < n_threads; ++j) {
				int i;
				worker_t *w = wo + j;
				for (i = 0; i < w->n; ++i)
					bits[w->buf[i]>>6] |= 1ull<<(w->buf[i]&0x3f);
				rest -= w->n;
			}
		}
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] %.3f million symbols remain; CPU time: %.3f; wall-clock / CPU: %.2f\n", __func__,
					rest / 1e6, cputime()-tcpu, (cputime()-tcpu)/(realtime()-treal));
	}
	free(wo); free(wr); free(tid);
	return bits;
}
