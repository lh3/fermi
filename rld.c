#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "rld.h"

#define RLD_IBITS_PLUS 4

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

rld_t *rld_init(int asize, int bbits)
{
	rld_t *e;
	e = calloc(1, sizeof(rld_t));
	e->n = 1;
	e->z = malloc(sizeof(void*));
	e->z[0] = calloc(RLD_LSIZE, 8);
	e->ssize = 1<<bbits;
	e->cnt = calloc(asize + 1, 8);
	e->mcnt = calloc(asize + 1, 8);
	e->abits = ilog2(asize) + 1;
	e->asize = asize;
	e->sbits = bbits;
	e->asize1 = asize + 1;
	e->offset0[0] = (e->asize1*16+63)/64;
	e->offset0[1] = (e->asize1*32+63)/64;
	return e;
}

void rld_destroy(rld_t *e)
{
	int i = 0;
	if (e == 0) return;
	for (i = 0; i < e->n; ++i) free(e->z[i]);
	free(e->z); free(e->cnt); free(e->mcnt); free(e->frame); free(e);
}

static inline void enc_next_block(rld_t *e, rlditr_t *itr)
{
	int i;
	if (itr->p + 1 - *itr->i == RLD_LSIZE) {
		++e->n;
		e->z = realloc(e->z, e->n * sizeof(void*));
		itr->i = e->z + e->n - 1;
		itr->shead = *itr->i = calloc(RLD_LSIZE, 8);
	} else itr->shead += e->ssize;
	if (e->cnt[0] - e->mcnt[0] >= 0x8000) {
		uint32_t *p = (uint32_t*)itr->shead;
		for (i = 0; i <= e->asize; ++i) p[i] = e->cnt[i] - e->mcnt[i];
		*p |= 1u<<31;
		itr->p = itr->shead + e->offset0[1];
	} else {
		uint16_t *p = (uint16_t*)itr->shead;
		for (i = 0; i <= e->asize; ++i) p[i] = e->cnt[i] - e->mcnt[i];
		itr->p = itr->shead + e->offset0[0];
	}
	itr->stail = itr->shead + e->ssize - 1;
	itr->r = 64;
	for (i = 0; i <= e->asize; ++i) e->mcnt[i] = e->cnt[i];
}

int rld_enc(rld_t *e, rlditr_t *itr, int l, uint8_t c)
{
	int w;
	uint64_t x = rld_delta_enc1(l, &w) << e->abits | c;
	w += e->abits;
	if (w >= itr->r && itr->p == itr->stail) enc_next_block(e, itr);
	if (w > itr->r) {
		w -= itr->r;
		*itr->p++ |= x >> w;
		*itr->p = x << (itr->r = 64 - w);
	} else itr->r -= w, *itr->p |= x << itr->r;
	e->cnt[0] += l;
	e->cnt[c + 1] += l;
	return 0;
}

void rld_rank_index(rld_t *e)
{
	uint64_t last, n_blks;

	n_blks = e->n_bytes * 8 / 64 / e->ssize + 1;
	last = rld_last_blk(e);
	{
		uint64_t i, k, *cnt;
		int j;
		cnt = alloca(e->asize * 8);
		e->ibits = ilog2(e->mcnt[0] / n_blks) + RLD_IBITS_PLUS;
		e->n_frames = ((e->mcnt[0] + (1<<e->ibits) - 1) >> e->ibits) + 1;
		e->frame = calloc(e->n_frames * e->asize1, 8);
		e->frame[0] = 0;
		for (j = 0; j < e->asize; ++j) cnt[j] = 0;
		for (i = e->ssize, k = 1; i <= last; i += e->ssize) {
			uint64_t sum, *p = rld_seek_blk(e, i);
			if (rld_size_bit(*p)) { // 32-bit count
				uint32_t *q = (uint32_t*)p;
				for (j = 1; j <= e->asize; ++j) cnt[j-1] += q[j];
			} else { // 16-bit count
				uint16_t *q = (uint16_t*)p;
				for (j = 1; j <= e->asize; ++j) cnt[j-1] += q[j];
			}
			for (j = 0, sum = 0; j < e->asize; ++j) sum += cnt[j];
			while (sum >= k<<e->ibits) ++k;
			if (k < e->n_frames) {
				uint64_t x = k * e->asize1;
				e->frame[x] = i;
				for (j = 0; j < e->asize; ++j) e->frame[x + j + 1] = cnt[j];
			}
		}
		assert(k >= e->n_frames - 1);
		for (k = 1; k < e->n_frames; ++k) { // fill zero cells
			uint64_t x = k * e->asize1;
			if (e->frame[x] == 0) {
				for (j = 0; j <= e->asize; ++j)
					e->frame[x + j] = e->frame[x - e->asize1 + j];
			}
		}
	}
}

uint64_t rld_enc_finish(rld_t *e, rlditr_t *itr)
{
	int i;
	enc_next_block(e, itr);
	e->n_bytes = (((uint64_t)(e->n - 1) * RLD_LSIZE) + (itr->p - *itr->i)) * 8;
	// recompute e->cnt as the accumulative count; e->mcnt[] keeps the marginal counts
	for (e->cnt[0] = 0, i = 1; i <= e->asize; ++i) e->cnt[i] += e->cnt[i - 1];
	rld_rank_index(e);
	return e->n_bytes;
}

int rld_dump(const rld_t *e, const char *fn)
{
	uint32_t a[2];
	uint64_t k;
	int i;
	FILE *fp;
	fp = strcmp(fn, "-")? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
	if (fp == 0) return -1;
	a[0] = e->asize; a[1] = e->sbits;
	fwrite("RLD\1", 1, 4, fp); // write magic
	fwrite(a, 4, 2, fp); // write asize and sbits
	fwrite(e->mcnt, 8, e->asize + 1, fp); // write the marginal counts
	fwrite(&e->n_bytes, 8, 1, fp);
	for (i = 0, k = e->n_bytes / 8; i < e->n - 1; ++i, k -= RLD_LSIZE)
		fwrite(e->z[i], 8, RLD_LSIZE, fp);
	fwrite(e->z[i], 8, k, fp);
	fwrite(&e->n_frames, 8, 1, fp);
	fwrite(e->frame, 8 * e->asize1, e->n_frames, fp); // FIXME: we cannot have >2GB frames with fwrite()
	fclose(fp);
	return 0;
}

rld_t *rld_restore(const char *fn)
{
	FILE *fp;
	rld_t *e;
	uint32_t a[2];
	char magic[5];
	uint64_t k, n_blks;
	int i;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	fread(magic, 1, 4, fp); magic[4] = 0;
	if (strcmp(magic, "RLD\1")) return 0; // read magic
	if (fread(a, 4, 2, fp) != 2) return 0; // read asize and sbits
	e = rld_init(a[0], a[1]);
	fread(e->mcnt, 8, e->asize + 1, fp);
	for (i = 0; i <= e->asize; ++i) e->cnt[i] = e->mcnt[i];
	for (e->cnt[0] = 0, i = 1; i <= e->asize; ++i) e->cnt[i] += e->cnt[i - 1];
	fread(&e->n_bytes, 8, 1, fp);
	if (e->n_bytes / 8 > RLD_LSIZE) { // allocate enough memory
		e->n = (e->n_bytes / 8 + RLD_LSIZE - 1) / RLD_LSIZE;
		e->z = realloc(e->z, e->n * sizeof(void*));
		for (i = 1; i < e->n; ++i)
			e->z[i] = calloc(RLD_LSIZE, 8);
	}
	for (i = 0, k = e->n_bytes / 8; i < e->n - 1; ++i, k -= RLD_LSIZE)
		fread(e->z[i], 8, RLD_LSIZE, fp);
	fread(e->z[i], 8, k, fp);
	fread(&e->n_frames, 8, 1, fp);
	e->frame = malloc(e->n_frames * e->asize1 * 8);
	fread(e->frame, 8 * e->asize1, e->n_frames, fp);
	fclose(fp);
	n_blks = e->n_bytes * 8 / 64 / e->ssize + 1;
	e->ibits = ilog2(e->mcnt[0] / n_blks) + RLD_IBITS_PLUS;
	return e;
}

uint64_t rld_rawlen(const rld_t *e)
{
	uint64_t x = rld_last_blk(e);
	return rld_seek_blk(e, x)[0];
}

#ifdef RLD_MAIN
int main(int argc, char *argv[])
{
	int i, l, c;
	rld_t *r = rld_enc_init(6, 3);
	for (i = 100000; i < 100100; ++i)
		rld_enc(r, i, 2);
	rld_enc_finish(r);
	rld_dec_init(r, 0);
	for (i = 100000; i < 100100; ++i) {
		l = rld_dec(r, &c);
		printf("%d, %d\n", l, c);
	}
	rld_destroy(r);
	return 0;
}
#endif
