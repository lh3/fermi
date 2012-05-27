#include <pthread.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "priv.h"

/**************************
 * Compute the bit vector *
 **************************/

#define BLOCK_SIZE 0x40000
#define TIMER_INTV 64

typedef struct {
	int start, step;
	const rld_t *e0, *e1;
	uint64_t *bits;
	int64_t *buf;
} worker_t;

static inline void update_bits(int n, const int64_t *buf, uint64_t *bits)
{
	const int64_t *q, *end = buf + n;
	for (q = buf; q != end; ++q) {
		uint64_t *p = bits + (*q>>6);
		uint64_t x = 1ull<<(*q&0x3f);
		__sync_or_and_fetch(p, x); // SEE ALSO: http://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Atomic-Builtins.html
	}
}

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int n = 0;
	int64_t i, k, x, n_processed = 0;
	uint64_t *ok;
	double tcpu, treal;
	tcpu = cputime(); treal = realtime();
	ok = alloca(8 * w->e0->asize);
	k = x = w->start;
	i = w->e0->mcnt[1] - 1;
	w->buf[n++] = i + k + 1;
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
			update_bits(n, w->buf, w->bits);
			if (fm_verbose >= 3 && ++n_processed % TIMER_INTV == 0)
				fprintf(stderr, "[M::%s@%d] processed %.3f million symbols in %.3f / %.1f seconds.\n", __func__, w->start,
						(double)n_processed*BLOCK_SIZE/1e6, cputime() - tcpu, (cputime() - tcpu) / (realtime() - treal));
			n = 0;
		}
		w->buf[n++] = k + i + 1;
	}
	if (n) update_bits(n, w->buf, w->bits);
	return 0;
}

uint64_t *fm_compute_gap_bits(const rld_t *e0, const rld_t *e1, int n_threads)
{
	uint64_t *bits;
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	int j;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e0 = e0; ww->e1 = e1;
		ww->step = n_threads;
		ww->start = j;
		ww->buf = xmalloc(BLOCK_SIZE * 8);
	}
	bits = xcalloc((e0->mcnt[0] + e1->mcnt[0] + 63) / 64, 8);
	for (j = 0; j < n_threads; ++j) w[j].bits = bits;
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	for (j = 0; j < n_threads; ++j) free(w[j].buf);
	free(w); free(tid);
	return bits;
}

/************************
 * Merge two FM-indexes *
 ************************/

rld_t *fm_merge(rld_t *e0, rld_t *e1, int n_threads)
{
	uint64_t i, n = e0->mcnt[0] + e1->mcnt[0], *bits;
	rlditr_t itr, itr0, itr1;
	rld_t *e;

	// compute the gap array
	bits = fm_compute_gap_bits(e0, e1, n_threads);
	free(e0->frame); free(e1->frame); // deallocate the rank indexes of e0 and e1; they are not needed any more
	e0->frame = e1->frame = 0;
	// initialize the FM-index to be returned, and all the three iterators
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr, 0);
	rld_itr_init(e0, &itr0, 0);
	rld_itr_init(e1, &itr1, 0);
	{
		uint64_t k = 1;
		int last = bits[0]&1;
		for (i = 1; i < n; ++i) {
			int c = bits[i>>6]>>(i&0x3f)&1;
			if (c != last) {
				if (last == 0) rld_dec_enc(e, &itr, e0, &itr0, k);
				else rld_dec_enc(e, &itr, e1, &itr1, k);
				last = c; k = 1;
			} else ++k;
		}
		if (k) {
			if (last == 0) rld_dec_enc(e, &itr, e0, &itr0, k);
			else rld_dec_enc(e, &itr, e1, &itr1, k);
		}
	}
	// finalize the merge
	assert(itr0.l == 0 && itr1.l == 0); // both e0 and e1 stream should be finished
	free(bits);
	rld_destroy(e0); rld_destroy(e1);
	rld_enc_finish(e, &itr);
	return e;
}

/**********************************
 * Append a string to an FM-index *
 **********************************/

rld_t *fm_merge_from_SA(rld_t *e0, int len, const uint8_t *T, const int *SA, const int64_t *rank_l)
{
	int64_t last = -1;
	int i;
	rlditr_t itr, itr0;
	rld_t *e;

	free(e0->frame); e0->frame = 0;
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr, 0);
	rld_itr_init(e0, &itr0, 0);
	for (i = 0; i < len; ++i) {
		if (rank_l[i] != last) {
			rld_dec_enc(e, &itr, e0, &itr0, rank_l[i] - last);
			last = rank_l[i];
		}
		rld_enc(e, &itr, 1, SA[i]? T[SA[i]-1] : 0);
	}
	if (last != e0->mcnt[0] - 1)
		rld_dec_enc(e, &itr, e0, &itr0, e0->mcnt[0] - 1 - last);
	rld_destroy(e0);
	rld_enc_finish(e, &itr);
	return e;
}

rld_t *fm_append(rld_t *e0, int len, const uint8_t *T)
{
	int c, k, *C, *p, *SA;
	uint64_t i, *oi, *rank_l;
	uint32_t *ws;
	rld_t *e;

	assert(T[len-1] == 0); // must be ended with a sentinel
	C = alloca(sizeof(int) * (e0->asize + 1));
	for (c = 0; c <= e0->asize; ++c) C[c] = 0;
	for (k = 0; k < len; ++k) ++C[T[k] + 1]; // marginal count
	for (c = 1; c <= e0->asize; ++c) C[c] += C[c-1]; // accumulative count
	ws = xmalloc((size_t)(len + 1) * 12);
	// set pointers
	rank_l = (uint64_t*)ws;
	SA = (int*)ws + 2 * (len + 1);
	// construct the suffix array
	ksa_sa(T, (int*)SA, len, e0->asize);
	// grab some memory from the stack
	p = alloca(sizeof(int) * e0->asize);
	oi = alloca(8 * e0->asize);
	// initialize the position array
	for (c = 0; c < e0->asize; ++c) p[c] = C[c+1] - 1; // point to the last element of each bucket
	// put the last sentinel 
	i = e0->mcnt[1] - 1;
	rank_l[p[0]--] = i;
	// compute the rank of long suffixes
	for (k = len - 2; k >= 0; --k) {
		if ((c = T[k]) != 0) {
			rld_rank1a(e0, i, oi);
			i = e0->cnt[c] + oi[c] - 1;
		} else i = e0->mcnt[1] - 1;
		rank_l[p[c]--] = i;
	}
	// sort the rank of long suffixes
	for (c = 1; c < e0->asize; ++c) // do not sort the sentinel bucket
		ks_introsort_uint64_t(C[c+1] - C[c], rank_l + C[c]);
	// merge to e0; e0 will be deallocated
	e = fm_merge_from_SA(e0, len, T, SA, (int64_t*)rank_l);
	free(ws);
	return e;
}
