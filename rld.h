#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define RLD_LBITS 23
#define RLD_LSIZE (1<<RLD_LBITS)
#define RLD_LMASK (RLD_LSIZE - 1)

typedef struct {
	int r; // bits remained in the last 64-bit integer
	uint64_t *p, *shead, *stail, **i;
	uint8_t *q;
} rlditr_t;

typedef struct __rld_t {
	// initialized in the constructor
	uint8_t asize, asize1; // alphabet size; asize1=asize+1
	int8_t abits; // bits required to store a symbol
	int8_t sbits; // bits per small block
	int8_t ibits; // modified during indexing; here for a better alignment
	int8_t offset0[2]; // 0 for 16-bit blocks; 1 for 32-bit blocks
	int ssize; // ssize = 1<<sbits
	// modified during encoding
	int n; // number of blocks (unchanged in decoding)
	uint64_t n_bytes; // total number of bits (unchanged in decoding)
	uint64_t **z; // the actual data (unchanged in decoding)
	uint64_t *cnt, *mcnt; // after enc_finish, cnt keeps the accumulative count and mcnt keeps the marginal
	// modified during indexing
	uint64_t n_frames;
	uint64_t *frame;
	// on-disk generation
	FILE *fp;
} rld_t;

#ifdef __cplusplus
extern "C" {
#endif

	rld_t *rld_init(int asize, int bbits, const char *fn);
	int rld_enc(rld_t *e, rlditr_t *itr, int64_t l, uint8_t c);
	uint64_t rld_enc_finish(rld_t *e, rlditr_t *itr);

	int rld_dump(const rld_t *e, const char *fn);
	rld_t *rld_restore(const char *fn);
	void rld_destroy(rld_t *e);

	uint64_t rld_rank11(const rld_t *e, uint64_t k, int c);
	void rld_rank1a(const rld_t *e, uint64_t k, uint64_t *ok);
	void rld_rank21(const rld_t *e, uint64_t k, uint64_t l, int c, uint64_t *ok, uint64_t *ol);
	void rld_rank2a(const rld_t *e, uint64_t k, uint64_t l, uint64_t *ok, uint64_t *ol);

#ifdef __cplusplus
}
#endif

#define rld_last_blk(e) ((e)->n_bytes>>3>>(e)->sbits<<(e)->sbits)
#define rld_seek_blk(e, k) ((e)->z[(k)>>RLD_LBITS] + ((k)&RLD_LMASK))

#define rld_size_bit(x) ((x)>>31&1)

static inline void rld_itr_init(const rld_t *e, rlditr_t *itr, uint64_t k)
{
	itr->i = e->z + (k >> RLD_LBITS);
	itr->shead = *itr->i + k%RLD_LSIZE;
	itr->stail = itr->shead + e->ssize - 1;
	itr->p = itr->shead + e->offset0[rld_size_bit(*itr->shead)];
	itr->q = (uint8_t*)itr->p;
	itr->r = 64;
}

#ifdef _USE_RLE6
static inline int64_t rld_dec0(const rld_t *r, rlditr_t *itr, int *c)
{
	int64_t l;
	*c = *itr->q >> 5;
	if (*c < 6) {
		l = *itr->q&0x1f;
		++itr->q;
	} else if (*c == 6) {
		*c = *itr->q>>2 & 0x7;
		if (*c != 6) {
			l = ((int64_t)itr->q[0]&0x3) | (int64_t)itr->q[1]<<2;
			itr->q += 2;
		} else {
			*c = (itr->q[0]&0x3) | itr->q[1]>>7<<2;
			l = ((uint32_t)itr->q[1]&0x7f) | (uint32_t)itr->q[2] << 7 | (uint32_t)itr->q[3] << 15;
			itr->q += 4;
		}
	} else return 0;
	return l;
}
#else
static inline int64_t rld_dec0(const rld_t *e, rlditr_t *itr, int *c)
{
	int w;
	uint64_t x;
	int64_t l, y = 0;
//	assert(itr->p <= itr->stail);
	x = itr->p[0] << (64 - itr->r) | (itr->p < itr->stail && itr->r < 64? itr->p[1] >> itr->r : 0);
	if (x>>63 == 0) {
		if ((w = 0x333333335555779bll>>(x>>59<<2)&0xf) == 0xb && x>>58 == 0) return 0;
		l = (x >> (64 - w)) - 1;
		y = x << w >> (64 - l) | 1u << l;
		w += l;
	} else w = y = 1;
	*c = x << w >> (64 - e->abits);
	w += e->abits;
	if (itr->r > w) itr->r -= w;
	else ++itr->p, itr->r = 64 + itr->r - w;
	return y;
}
#endif

static inline int64_t rld_dec(const rld_t *e, rlditr_t *itr, int *_c)
{
	int64_t l = rld_dec0(e, itr, _c);
	if (l == 0) {
		uint64_t last = rld_last_blk(e);
		if (itr->p - *itr->i > RLD_LSIZE - e->ssize) itr->shead = *++itr->i;
		else itr->shead += e->ssize;
		if (itr->shead == rld_seek_blk(e, last)) return -1;
		itr->p = itr->shead + e->offset0[rld_size_bit(*itr->shead)];
		itr->q = (uint8_t*)itr->p;
		itr->stail = itr->shead + e->ssize - 1;
		itr->r = 64;
		return rld_dec0(e, itr, _c);
	} else return l;
}

#endif
