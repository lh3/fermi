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
	e->ssize = 1<<bbits;
	e->stail = e->shead + e->ssize - 1;
	e->cnt = calloc(asize + 1, 8);
	e->mcnt = calloc(asize + 1, 8);
	e->abits = ilog2(asize) + 1;
	e->asize = asize;
	e->sbits = bbits;
	e->asize1 = asize + 1;
	e->o0[0] = (e->asize1*16+63)/64;
	e->o0[1] = (e->asize1*32+63)/64;
	e->r = 64;
	e->p = e->shead + e->o0[0];
	return e;
}

static inline void enc_next_block(rld_t *e)
{
	int i;
	if (e->p + 1 - e->lhead == RLD_LSIZE) {
		++e->n;
		e->z = realloc(e->z, e->n * sizeof(void*));
		e->lhead = e->shead = e->z[e->n - 1] = calloc(RLD_LSIZE, 8);
	} else e->shead += e->ssize;
	if (e->cnt[0] - e->mcnt[0] >= 0x8000) {
		uint32_t *p = (uint32_t*)e->shead;
		for (i = 0; i <= e->asize; ++i) p[i] = e->cnt[i] - e->mcnt[i];
		*p |= 1u<<31;
		e->p = e->shead + e->o0[1];
	} else {
		uint16_t *p = (uint16_t*)e->shead;
		for (i = 0; i <= e->asize; ++i) p[i] = e->cnt[i] - e->mcnt[i];
		e->p = e->shead + e->o0[0];
	}
	e->stail = e->shead + e->ssize - 1;
	e->r = 64;
	for (i = 0; i <= e->asize; ++i) e->mcnt[i] = e->cnt[i];
}

int rld_enc(rld_t *e, int l, uint8_t c)
{
	int w;
	uint64_t x = rld_delta_enc1(l, &w) << e->abits | c;
	w += e->abits;
	if (w > e->r && e->p == e->stail) enc_next_block(e);
	if (w > e->r) {
		w -= e->r;
		*e->p++ |= x >> w;
		*e->p = x << (e->r = 64 - w);
	} else e->r -= w, *e->p |= x << e->r;
	e->cnt[0] += l;
	e->cnt[c + 1] += l;
	return 0;
}

uint64_t rld_enc_finish(rld_t *e)
{
	int i;
	uint64_t last;
	enc_next_block(e);
	e->n_bits = (((uint64_t)(e->n - 1) * RLD_LSIZE) + (e->p - e->lhead)) * 64;
	// recompute e->cnt as the accumulative count; e->mcnt[] keeps the marginal counts
	for (e->cnt[0] = 0, i = 1; i <= e->asize; ++i) e->cnt[i] += e->cnt[i - 1];
	last = rld_last_blk(e);
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
	uint64_t last, n_blks;

	idx = calloc(1, sizeof(rldidx_t));
	n_blks = e->n_bits / 64 / e->ssize + 1;
	last = rld_last_blk(e);
	{
		uint64_t i, k, *cnt;
		int j;
		cnt = alloca(e->asize * 8);
		idx->b = ilog2(e->mcnt[0] / n_blks) + 3;
		idx->n = ((e->mcnt[0] + (1<<idx->b) - 1) >> idx->b) + 1;
		idx->s = calloc(idx->n * e->asize1, 8);
		idx->s[0] = 0;
		for (j = 0; j < e->asize; ++j) cnt[j] = 0;
		for (i = e->ssize, k = 1; i <= last; i += e->ssize) {
			uint64_t sum, *p = rld_seek_blk(e, i);
			if (*p>>63) { // 32-bit count
				uint32_t *q = (uint32_t*)p;
				for (j = 1; j <= e->asize; ++j) cnt[j-1] += q[j];
			} else { // 16-bit count
				uint16_t *q = (uint16_t*)p;
				for (j = 1; j <= e->asize; ++j) cnt[j-1] += q[j];
			}
			for (j = 0, sum = 0; j < e->asize; ++j) sum += cnt[j];
			while (sum >= k<<idx->b) ++k;
			if (k < idx->n) {
				uint64_t x = k * e->asize1;
				idx->s[x] = i;
				for (j = 0; j < e->asize; ++j) idx->s[x + j - 1] = cnt[j];
			}
		}
		for (j = 0; j < e->asize; ++j) printf("[0] %d, %lld, %lld\n", j, cnt[j], e->mcnt[j+1]);
		assert(k >= idx->n - 1);
		for (k = 1; k < idx->n; ++k) { // fill zero cells
			uint64_t x = k * e->asize1;
			if (idx->s[x] == 0) {
				for (j = 0; j <= e->asize; ++j)
					idx->s[x + j] = idx->s[x - e->asize1 + j];
			}
		}
	}
	return idx;
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
