#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>

#define RLD_LBITS 23
#define RLD_LSIZE (1<<RLD_LBITS)
#define RLD_LMASK (RLD_LSIZE - 1)

typedef struct {
	// static members (initialized in the constructor)
	int asize; // alphabet size
	int abits; // bits required to store a symbol
	int sbits; // bits per small block
	int ssize; // ssize = 1<<sbits
	// dynamic members
	int n; // number of blocks (unchanged in decoding)
	int r; // bits remaining in the last 64-bit integer
	uint64_t n_bits; // total number of bits (unchanged in decoding)
	uint64_t **z; // the actual data (unchanged in decoding)
	uint64_t *cnt;
	uint64_t *p, *lhead, *shead, *stail;
} rld_t;

typedef struct {
	int rsize, ibits;
	uint64_t *r;
} rldidx_t;

#ifdef __cplusplus
extern "C" {
#endif

	rld_t *rld_enc_init(int asize, int bbits);
	int rld_enc(rld_t *e, int l, uint8_t c);
	uint64_t rld_enc_finish(rld_t *e);
	uint64_t rld_rawlen(const rld_t *e);

#ifdef __cplusplus
}
#endif

#define rld_last_blk(e) ((e)->n_bits>>6>>(e)->sbits<<(e)->sbits)
#define rld_seek_blk(e, k) ((e)->z[(k)>>RLD_LBITS] + ((k)&RLD_LMASK))

static inline void rld_dec_init(rld_t *e, uint64_t k)
{
	e->lhead = e->z[k>>RLD_LBITS];
	e->shead = e->lhead + k%RLD_LSIZE;
	e->stail = e->shead + e->ssize - 1;
	e->p = e->shead + e->asize;
	e->r = 64;
}

static inline uint32_t rld_dec0(rld_t *e)
{
	int y = 0, w, l;
	uint64_t x;
	x = e->p[0] << (64 - e->r) | (e->p < e->stail && e->r < 64? e->p[1] >> e->r : 0);
	if (x>>63 == 0) {
		if ((w = 0x333333335555779bll>>(x>>59<<2)&0xf) == 0xb && x>>58 == 0) return 0;
		l = (x >> (64 - w)) - 1;
		y = x << w >> (64 - l) | 1u << l;
		w += l;
	} else w = y = 1;
	y = y << e->abits | x << w >> (64 - e->abits);
	w += e->abits;
	if (e->r > w) e->r -= w;
	else ++e->p, e->r = 64 + e->r - w;
	return y;
}

static inline uint32_t rld_dec(rld_t *e)
{
	uint32_t c = rld_dec0(e);
	if (c == 0) {
		if (e->p - e->lhead > RLD_LSIZE - e->ssize) e->shead = ++e->lhead;
		else e->shead += e->ssize;
		e->p = e->shead + e->asize;
		e->stail = e->shead + e->ssize - 1;
		e->r = 64;
		return rld_dec0(e);
	} else return c;
}
/*
static inline uint64_t rld_rank1(rld_t *e, const rldidx_t *r, uint64_t k, int c)
{
	uint64_t x = r->r[k>>r->ibits], *end, y;
	uint32_t a;
	int i;
	e->p = rld_seek_blk(e, x); e->r = 64;
	end = (x + RLD_LMASK) >> RLD_LBITS << RLD_LBITS;
	while (*e->p > k) {
		if ((x += e->ssize) == end) e->p = ld_seek_blk(e, end);
		else e->p += e->ssize;
	}
	x = *e->p;
	if (c == 0) {
		int i;
		for (i = 1, y = *e->p; i < e->asize; ++i) y -= e->p[i];
	} else y = e->p[c];
	return y;
}
*/
#endif
