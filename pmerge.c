#include <pthread.h>
#include <stdio.h>
#include "fermi.h"
#include "rld.h"
#include "utils.h"

/* This file gives two implementations of computing the gap array for the BWT
 * merge. The first implementation requires a processor-specific instruction
 * for the or-and-fetch atomic memory access [1], which is built in in recent
 * versions of GCC as __sync_or_and_fetch() [2]. With this built-in function we
 * can achieve multithreading in a lock-free manner and simplify coding as
 * well.
 *
 * For processors do not have such an atomic instruction or for compilers that
 * do not have the __sync_or_and_fetch() builtin, we can use the second
 * implementation, which in its current form implicitly sets a barrier before
 * we start to write the critical section. This implementation is less
 * efficient but more portable.
 *
 * [1] http://en.wikipedia.org/wiki/Fetch-and-add
 * [2] http://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Atomic-Builtins.html
 */

#ifndef _NO_ATOMIC

#define BLOCK_SIZE 0x40000
#define TIMER_INTV 64

typedef struct {
	int start, step;
	const rld_t *e0, *e1;
	uint64_t *bits;
} worker_t;

static inline void update_bits(int n, const int64_t *buf, uint64_t *bits)
{
	const int64_t *q, *end = buf + n;
	for (q = buf; q != end; ++q) {
		uint64_t *p = bits + (*q>>6);
		uint64_t x = 1ull<<(*q&0x3f);
		__sync_or_and_fetch(p, x); // IMPORTANT: *p|=x is not atomic
	}
}

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int n = 0;
	int64_t i, k, x, *buf, n_processed = 0;
	uint64_t *ok;
	double tcpu, treal;
	tcpu = cputime(); treal = realtime();
	ok = alloca(8 * w->e0->asize);
	buf = malloc(BLOCK_SIZE * 8);
	k = x = w->start;
	i = w->e0->mcnt[1] - 1;
	buf[n++] = i + k + 1;
	for (;;) {
		int c = rld_rank1a(w->e1, k, ok);
		if (c == 0) {
			x += w->step;
			if (x >= w->e1->mcnt[1]) break;
			k = x;
			i = w->e0->mcnt[1] - 1;
		} else {
			k = w->e1->cnt[c] + ok[c] - 1;
			rld_rank1a(w->e0, i, ok);
			i = w->e0->cnt[c] + ok[c] - 1;
		}
		if (n == BLOCK_SIZE) {
			update_bits(n, buf, w->bits);
			if (fm_verbose >= 3 && ++n_processed % TIMER_INTV == 0)
				fprintf(stderr, "[M::%s@%d] processed %f million symbols in %.3f * %.1f seconds.\n", __func__, w->start,
						(double)n_processed*BLOCK_SIZE/1e6, cputime() - tcpu, (cputime() - tcpu) / (realtime() - treal));
			n = 0;
		}
		buf[n++] = k + i + 1;
	}
	if (n) update_bits(n, buf, w->bits);
	free(buf);
	return 0;
}

uint64_t *fm_compute_gap_bits(const rld_t *e0, const rld_t *e1, int n_threads)
{
	uint64_t *bits;
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	int j;

	bits = xcalloc((e0->mcnt[0] + e1->mcnt[0] + 63) / 64, 8);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e0 = e0; ww->e1 = e1;
		ww->step = n_threads;
		ww->start = j;
		ww->bits = bits;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(w); free(tid);
	return bits;
}

#else

#define BLOCK_SIZE 0x1000000
#define PWRITE_MIN 4

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef struct {
	int step, n;
	int64_t size, i, k, x, *buf;
	const rld_t *e0, *e1;
} worker_t;

typedef struct {
	int tid, n_workers, mask;
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
	for (j = 0; j < w->n_workers; ++j) {
		worker_t *wo = w->w + j;
		for (i = 0; i < wo->n; ++i) {
			uint64_t x = wo->buf[i];
			if ((x>>6&w->mask) == w->tid)
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
	int j, n_writers;
	double tcpu, treal;

	bits = xcalloc((e0->mcnt[0] + e1->mcnt[0] + 63) / 64, 8);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	wo = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (j = 0; j < n_threads; ++j) {
		worker_t *w = wo + j;
		w->e0 = e0; w->e1 = e1;
		w->step = n_threads;
		w->size = ((BLOCK_SIZE < rest? BLOCK_SIZE : rest) + n_threads - 1) / n_threads;
		w->k = w->x = j;
		w->i = w->e0->mcnt[1] - 1;
		w->buf = xmalloc(w->size * 8);
	}
	n_writers = n_threads;
	kroundup32(n_writers);
	if (n_writers > n_threads) n_writers >>= 1;
	wr = (writer_t*)calloc(n_writers, sizeof(writer_t));
	for (j = 0; j < n_writers; ++j) {
		writer_t *r = wr + j;
		r->n_workers = n_threads;
		r->mask = n_writers - 1;
		r->tid = j;
		r->bits = bits;
		r->w = wo;
	}
	tcpu = cputime(); treal = realtime();
	while (rest) {
		for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, wo + j);
		for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
		if (n_threads >= PWRITE_MIN) {
			for (j = 0; j < n_writers; ++j) pthread_create(&tid[j], &attr, writer, wr + j);
			for (j = 0; j < n_writers; ++j) pthread_join(tid[j], 0);
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
	for (j = 0; j < n_threads; ++j) free(wo[j].buf);
	free(wo); free(wr); free(tid);
	return bits;
}
#endif
