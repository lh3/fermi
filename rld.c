#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "rld.h"

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2(uint32_t v)
{
	register uint32_t t, tt;
	if ((tt = (v >> 16))) return (t = (tt >> 8)) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = (v >> 8)) ? 8 + LogTable256[t] : LogTable256[v];
}

inline uint32_t rld_delta_enc1(uint32_t x, int *width)
{
	int y = ilog2(x);
	int z = ilog2(y + 1);
	*width = (z<<1) + 1 + y;
	return (x^(1<<y)) | (y+1)<<y;
}

rld_t *rld_enc_init(int asize, int bbits)
{
	rld_t *e;
	e = calloc(1, sizeof(rld_t));
	e->z = malloc(sizeof(void*));
	e->z[0] = calloc(RLD_LSIZE, 8);
	e->n = 1;
	e->shead = e->lhead = e->z[0];
	e->p = e->lhead + asize + 1;
	e->ssize = 1<<bbits;
	e->stail = e->shead + e->ssize - 1;
	e->cnt = calloc(asize + 1, 8);
	e->abits = ilog2(asize) + 1;
	e->r = 64;
	e->asize = asize;
	e->sbits = bbits;
	return e;
}

static inline void enc_next_block(rld_t *e, int r)
{
	int i;
	if (e->p + 1 - e->lhead == RLD_LSIZE) {
		++e->n;
		e->z = realloc(e->z, e->n * sizeof(void*));
		e->lhead = e->shead = e->z[e->n - 1] = calloc(RLD_LSIZE, 8);
	} else e->shead += e->ssize;
	for (i = 0; i <= e->asize; ++i) e->shead[i] = e->cnt[i];
	e->stail = e->shead + e->ssize - 1;
	e->p = e->shead + e->asize + 1;
	e->r = r;
}

int rld_enc(rld_t *e, int l, uint8_t c)
{
	int w;
	uint64_t x = rld_delta_enc1(l, &w) << e->abits | c;
	w += e->abits;
	if (w > e->r) {
		if (e->p == e->stail) { // jump to the next block
			enc_next_block(e, 64 - w);
			*e->p |= x << e->r;
		} else {
			w -= e->r;
			*e->p++ |= x >> w;
			*e->p = x << (e->r = 64 - w);
		}
	} else e->r -= w, *e->p |= x << e->r;
	e->cnt[0] += l;
	e->cnt[c + 1] += l;
	return 0;
}

uint64_t rld_enc_finish(rld_t *e)
{
	int i;
	enc_next_block(e, 64);
	e->n_bits = (((uint64_t)(e->n - 1) * RLD_LSIZE) + (e->p - e->lhead)) * 64;
	// recompute e->cnt as the accumulative count
	for (e->cnt[0] = 0, i = 1; i <= e->asize; ++i) e->cnt[i] += e->cnt[i - 1];
	return e->n_bits;
}

uint64_t rld_rawlen(const rld_t *e)
{
	uint64_t x = rld_last_blk(e);
	return rld_seek_blk(e, x)[0];
}

rldidx_t *rld_index(const rld_t *e)
{
	rldidx_t *idx;
	uint64_t i, last, k, rawlen, n_blks;
	idx = calloc(1, sizeof(rldidx_t));
	n_blks = e->n_bits / 64 / e->ssize + 1;
	last = rld_last_blk(e);
	rawlen = rld_seek_blk(e, last)[0];
	idx->ibits = ilog2(rawlen / n_blks) + 1;
	idx->rsize = ((rawlen + (1<<idx->ibits) - 1) >> idx->ibits) + 1;
	idx->r = calloc(idx->rsize, 8);
	idx->r[0] = 0;
	for (i = e->ssize, k = 1; i <= last; i += e->ssize) {
		uint64_t x = *rld_seek_blk(e, i);
		while (x >= k<<idx->ibits) ++k;
		if (k < idx->rsize) idx->r[k] = i;
	}
	assert(k >= idx->rsize - 1);
	for (k = 1, last = 0; k < idx->rsize; ++k) {
		if (idx->r[k]) last = idx->r[k];
		else idx->r[k] = last;
	}
	//for (k = 0; k < idx->rsize; ++k)
	//	if (*rld_seek_blk(e, idx->r[k]) > k<<idx->ibits)
	//		printf("%lld > %lld\n", *rld_seek_blk(e, idx->r[k]), k<<idx->ibits);
	return idx;
}

uint64_t rld_backward_search(rld_t *e, const rldidx_t *r, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end)
{
	uint64_t k, l, ok, ol;
	int i, c;
	c = str[len - 1];
	k = e->cnt[c]; l = e->cnt[c + 1] - 1;
	for (i = len - 2; i >= 0; --i) {
		c = str[i];
		rld_rank12(e, r, k - 1, l, c, &ok, &ol);
		k = e->cnt[c] + ok;
		l = e->cnt[c] + ol - 1;
		if (k > l) break;
	}
	if (k > l) return 0;
	*sa_beg = k; *sa_end = l;
	return l - k + 1;
}

#ifdef RLD_MAIN
int main(int argc, char *argv[])
{
	int k, i, n = 100000, N = 500, a = 10, b = 1;
	rld_t *r = rld_enc_init(6, 5);
	for (i = 1; i < n; ++i)
		rld_enc(r, i%a+b, 0);
	fprintf(stderr, "# bytes: %f\n", rld_enc_finish(r) / 8.);
	for (k = 0; k < N; ++k) {
		int j = 0;
		rld_dec_init(r, j);
		for (i = 1; i < n; ++i)
			if (i%a+b != rld_dec(r) >> 3)
				fprintf(stderr, "Bug!\n");
	}
	for (i = 1; i < r->n; ++i) free(r->z);
	free(r->cnt); free(r->z); free(r);
	return 0;
}
#endif
