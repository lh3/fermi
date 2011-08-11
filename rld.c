#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "rld.h"

#define RLD_IBITS_PLUS 3

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

static void dump_header(const rld_t *e, FILE *fp)
{
	int32_t a[2];
	a[0] = RLD_ASIZE; a[1] = e->sbits;
	fwrite("RLD\1", 1, 4, fp); // write magic
	fwrite(a, 4, 2, fp); // write asize and sbits
	fwrite(e->mcnt, 8, RLD_ASIZE + 1, fp); // write the marginal counts
	fwrite(&e->n_bytes, 8, 1, fp);
}

rld_t *rld_init(int asize, int bbits, const char *fn)
{
	rld_t *e;
	e = calloc(1, sizeof(rld_t));
	e->n = 1;
	e->z = malloc(sizeof(void*));
	e->z[0] = calloc(RLD_LSIZE, 8);
	e->ssize = 1<<bbits;
	e->cnt = calloc(RLD_ASIZE + 1, 8);
	e->mcnt = calloc(RLD_ASIZE + 1, 8);
	e->sbits = bbits;
	if (fn) {
		e->fp = strcmp(fn, "-")? fopen(fn, "wb+") : stdout;
		dump_header(e, e->fp);
	}
	return e;
}

void rld_destroy(rld_t *e)
{
	int i = 0;
	if (e == 0) return;
	for (i = 0; i < e->n; ++i) free(e->z[i]);
	if (e->fp) fclose(e->fp);
	free(e->z); free(e->cnt); free(e->mcnt); free(e->frame); free(e);
}

static inline void enc_next_block(rld_t *e, rlditr_t *itr)
{
	int i;
	uint16_t *p;
	if (itr->q + 1 - *itr->i == RLD_LSIZE) {
		if (e->fp) {
			fwrite(e->z[e->n - 1], 1, RLD_LSIZE, e->fp);
			free(e->z[e->n - 1]);
			e->z[e->n - 1] = 0;
		}
		++e->n;
		e->z = realloc(e->z, e->n * sizeof(void*));
		itr->i = e->z + e->n - 1;
		itr->shead = *itr->i = calloc(RLD_LSIZE, 1);
	} else itr->shead += e->ssize;
	for (p = (uint16_t*)itr->shead, i = 0; i <= RLD_ASIZE; ++i)
		p[i] = e->cnt[i] - e->mcnt[i];
	itr->q = itr->shead + RLD_OFFSET;
	itr->stail = itr->shead + e->ssize - 1;
	for (i = 0; i <= RLD_ASIZE; ++i) e->mcnt[i] = e->cnt[i];
}

inline void rld_enc0(rld_t *r, rlditr_t *itr, int l, uint8_t c)
{
	if (itr->q + 1 > (uint8_t*)itr->stail) {
		*itr->q = 0xff;
		enc_next_block(r, itr);
	}
	*itr->q++ = c<<5 | l;
	r->cnt[0] += l;
	r->cnt[c + 1] += l;
}

int rld_enc(rld_t *r, rlditr_t *itr, int l, uint8_t c)
{
	const int max_l = (1<<5) - 1;
	for (; l > max_l; l -= max_l)
		rld_enc0(r, itr, max_l, c);
	rld_enc0(r, itr, l, c);
	return 0;
}

void rld_rank_index(rld_t *e)
{
	uint64_t last, n_blks;

	n_blks = e->n_bytes / e->ssize + 1;
	last = rld_last_blk(e);
	{
		uint64_t i, k, *cnt;
		int j;
		cnt = alloca(RLD_ASIZE * 8);
		e->ibits = ilog2(e->mcnt[0] / n_blks) + RLD_IBITS_PLUS;
		e->n_frames = ((e->mcnt[0] + (1<<e->ibits) - 1) >> e->ibits) + 1;
		e->frame = calloc(e->n_frames * RLD_ASIZE1, 8);
		e->frame[0] = 0;
		for (j = 0; j < RLD_ASIZE; ++j) cnt[j] = 0;
		for (i = e->ssize, k = 1; i <= last; i += e->ssize) {
			uint64_t sum;
			uint16_t *p = (uint16_t*)rld_seek_blk(e, i);
			for (j = 1; j <= RLD_ASIZE; ++j) cnt[j-1] += p[j];
			for (j = 0, sum = 0; j < RLD_ASIZE; ++j) sum += cnt[j];
			while (sum >= k<<e->ibits) ++k;
			if (k < e->n_frames) {
				uint64_t x = k * RLD_ASIZE1;
				e->frame[x] = i;
				for (j = 0; j < RLD_ASIZE; ++j) e->frame[x + j + 1] = cnt[j];
			}
		}
		assert(k >= e->n_frames - 1);
		for (k = 1; k < e->n_frames; ++k) { // fill zero cells
			uint64_t x = k * RLD_ASIZE1;
			if (e->frame[x] == 0) {
				for (j = 0; j <= RLD_ASIZE; ++j)
					e->frame[x + j] = e->frame[x - RLD_ASIZE1 + j];
			}
		}
	}
}

uint64_t rld_enc_finish(rld_t *e, rlditr_t *itr)
{
	int i;
	enc_next_block(e, itr);
	e->n_bytes = (((uint64_t)(e->n - 1) * RLD_LSIZE) + (itr->q - *itr->i));
	// recompute e->cnt as the accumulative count; e->mcnt[] keeps the marginal counts
	for (e->cnt[0] = 0, i = 1; i <= RLD_ASIZE; ++i) e->cnt[i] += e->cnt[i - 1];
	if (e->fp) {
		uint64_t k;
		fseek(e->fp, 28 + RLD_ASIZE, SEEK_SET);
		for (i = 0, k = e->n_bytes; i < e->n - 1; ++i, k -= RLD_LSIZE) {
			e->z[i] = malloc(RLD_LSIZE);
			fread(e->z[i], 1, RLD_LSIZE, e->fp);
		}
		fwrite(e->z[i], 1, k, e->fp);
		rld_rank_index(e);
		fwrite(&e->n_frames, 8, 1, e->fp);
		fwrite(e->frame, 8 * RLD_ASIZE1, e->n_frames, e->fp);
		rewind(e->fp);
		dump_header(e, e->fp);
		fclose(e->fp);
		e->fp = 0;
	} else rld_rank_index(e);
	return e->n_bytes;
}

int rld_dump(const rld_t *e, const char *fn)
{
	uint64_t k;
	int i;
	FILE *fp;
	fp = strcmp(fn, "-")? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
	if (fp == 0) return -1;
	dump_header(e, fp);
	for (i = 0, k = e->n_bytes; i < e->n - 1; ++i, k -= RLD_LSIZE)
		fwrite(e->z[i], 1, RLD_LSIZE, fp);
	fwrite(e->z[i], 1, k, fp);
	fwrite(&e->n_frames, 8, 1, fp);
	fwrite(e->frame, 8 * RLD_ASIZE1, e->n_frames, fp);
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
	e = rld_init(a[0], a[1], 0);
	fread(e->mcnt, 8, RLD_ASIZE + 1, fp);
	for (i = 0; i <= RLD_ASIZE; ++i) e->cnt[i] = e->mcnt[i];
	for (e->cnt[0] = 0, i = 1; i <= RLD_ASIZE; ++i) e->cnt[i] += e->cnt[i - 1];
	fread(&e->n_bytes, 8, 1, fp);
	if (e->n_bytes > RLD_LSIZE) { // allocate enough memory
		e->n = (e->n_bytes + RLD_LSIZE - 1) / RLD_LSIZE;
		e->z = realloc(e->z, e->n * sizeof(void*));
		for (i = 1; i < e->n; ++i)
			e->z[i] = calloc(RLD_LSIZE, 1);
	}
	for (i = 0, k = e->n_bytes; i < e->n - 1; ++i, k -= RLD_LSIZE)
		fread(e->z[i], 1, RLD_LSIZE, fp);
	fread(e->z[i], 1, k, fp);
	fread(&e->n_frames, 8, 1, fp);
	e->frame = malloc(e->n_frames * RLD_ASIZE1 * 8);
	fread(e->frame, 8 * RLD_ASIZE1, e->n_frames, fp);
	fclose(fp);
	n_blks = e->n_bytes / e->ssize + 1;
	e->ibits = ilog2(e->mcnt[0] / n_blks) + RLD_IBITS_PLUS;
	return e;
}

uint64_t rld_rawlen(const rld_t *e)
{
	uint64_t x = rld_last_blk(e);
	return rld_seek_blk(e, x)[0];
}

static inline uint64_t rld_locate_blk(const rld_t *e, rlditr_t *itr, uint64_t k, uint64_t *cnt, uint64_t *sum)
{
	int j;
	uint64_t c = 0, *z = e->frame + (k>>e->ibits) * RLD_ASIZE1;
	uint16_t *p;
	uint8_t *q;
	itr->i = e->z + (*z>>RLD_LBITS);
	q = itr->q = *itr->i + (*z&RLD_LMASK);
	for (j = 1, *sum = 0; j < RLD_ASIZE1; ++j) *sum += (cnt[j-1] = z[j]);
	while (1) { // seek to the small block
		q += e->ssize;
		if (q - *itr->i == RLD_LSIZE) q = *++itr->i;
		c = *(uint16_t*)q;
		if (*sum + c > k) break;
		for (p = (uint16_t*)q, j = 0; j < RLD_ASIZE; ++j) cnt[j] += p[j+1];
		*sum += c;
		itr->q = q;
	}
	itr->shead = itr->q;
	itr->stail = itr->shead + e->ssize - 1;
	itr->q += RLD_OFFSET;
	return c + *sum;
}

uint64_t rld_rank11(const rld_t *e, uint64_t k, int c)
{
	uint64_t y, z, *cnt;
	rlditr_t itr;
	if (k == (uint64_t)-1) return 0;
	cnt = alloca(RLD_ASIZE1 * 8);
	rld_locate_blk(e, &itr, k, cnt, &z);
	y = cnt[c];
	++k; // because k is the coordinate but not length
	while (1) {
		int a = -1, l;
		l = rld_dec0(e, &itr, &a);
		if (z + l >= k) return y + (a == c? k - z: 0);
		z += l;
		if (a == c) y += l;
	}
}

void rld_rank21(const rld_t *e, uint64_t k, uint64_t l, int c, uint64_t *ok, uint64_t *ol) // FIXME: can be faster
{
	*ok = rld_rank11(e, k, c);
	*ol = rld_rank11(e, l, c);
}

void rld_rank1a(const rld_t *e, uint64_t k, uint64_t *ok)
{
	uint64_t z, l;
	int a = -1;
	rlditr_t itr;
	if (k == (uint64_t)-1) {
		int a;
		for (a = 0; a < RLD_ASIZE; ++a) ok[a] = 0;
		return;
	}
	rld_locate_blk(e, &itr, k, ok, &z);
	++k; // because k is the coordinate but not length
	while (1) {
		l = rld_dec0(e, &itr, &a);
		if (z + l >= k) break;
		z += l; ok[a] += l;
	}
	ok[a] += k - z;
}

void rld_rank2a(const rld_t *e, uint64_t k, uint64_t l, uint64_t *ok, uint64_t *ol)
{
	uint64_t z, y, len;
	rlditr_t itr;
	int a = -1;
	if (k == (uint64_t)-1) { // special treatment for k==-1
		for (a = 0; a < RLD_ASIZE; ++a) ok[a] = 0;
		rld_rank1a(e, l, ol);
		return;
	}
	y = rld_locate_blk(e, &itr, k, ok, &z); // locate the block bracketing k
	++k; // because k is the coordinate but not length
	while (1) { // compute ok[]
		len = rld_dec0(e, &itr, &a);
		if (z + len >= k) break;
		z += len; ok[a] += len;
	}
	if (y > l) { // we do not need to decode other blocks
		int b;
		++l; // for a similar reason to ++l
		for (b = 0; b < RLD_ASIZE; ++b) ol[b] = ok[b]; // copy ok[] to ol[]
		ok[a] += k - z; // finalize ok[a]
		if (z + len < l) { // we need to decode the next run
			z += len; ol[a] += len;
			while (1) {
				len = rld_dec0(e, &itr, &a);
				if (z + len >= l) break;
				z += len; ol[a] += len;
			}
		}
		ol[a] += l - z;
	} else { // we have to decode other blocks
		ok[a] += k - z;
		rld_rank1a(e, l, ol);
	}
}

#ifdef RLD_MAIN
int main(int argc, char *argv[])
{
	int i, l, c;
	rlditr_t itr;
	rld_t *r = rld_init(6, 6, 0);
	rld_itr_init(r, &itr, 0);
	for (i = 1; i < 63; ++i)
		rld_enc(r, &itr, i, 2);
	rld_enc_finish(r, &itr);
	rld_itr_init(r, &itr, 0);
	for (i = 1; i < 63; ++i) {
		l = rld_dec(r, &itr, &c);
		printf("%d, %d\n", l, c);
	}
	rld_destroy(r);
	return 0;
}
#endif
