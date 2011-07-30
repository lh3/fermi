#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>

#define RLD_LBITS 23
#define RLD_LSIZE (1<<RLD_LBITS)

#define rld_blk_size(bbits, balpha) (((1<<(bbits)) << 6) / ((balpha) + 1))

typedef struct {
	// static members
	int asize; // alphabet size
	int abits; // bits required to store a symbol
	int sbits; // bits per small block
	int ssize; // ssize = 1<<sbits
	// dynamic members in encoding
	int n; // number of blocks
	// dynamic members in decoding
	int r; // bits remaining in the last 64-bit integer
	uint64_t *cnt;
	uint64_t **z;
	uint64_t *p, *lhead, *shead, *stail;
} rld_t;

#ifdef __cplusplus
extern "C" {
#endif

	rld_t *rld_enc_init(int asize, int bbits);
	int rld_enc(rld_t *e, int l, uint8_t c);
	uint64_t rld_enc_finish(rld_t *e);

#ifdef __cplusplus
}
#endif

extern int16_t rld_ddec_table[256];

static inline int rld_delta_dec1(uint64_t x, int *w)
{
	if (x >> 63) {
		*w = 1;
		return 1;
	} else {
		int z = rld_ddec_table[x >> 55];
		if (z > 0) {
			int a = z>>4, b = z&0xf;
			*w = a + b;
			return x << b >> (64 - a) | 1u << a;
		} else if (z < 0) {
			z = ~z;
			*w = z & 0xf;
			return z>>4;
		}
	}
	return 0;
}

static inline void rld_dec_init(rld_t *e, uint64_t k)
{
	e->lhead = e->z[k>>RLD_LBITS];
	e->shead = e->lhead + k%RLD_LSIZE;
	e->stail = e->shead + e->ssize - 1;
	e->p = e->shead + e->asize;
	e->r = 64;
}

static inline uint32_t rld_dec0a(rld_t *e)
{
	int y = 0, w;
	uint64_t x;
	x = e->p[0] << (64 - e->r) | (e->p < e->stail && e->r < 64? e->p[1] >> e->r : 0);
	/* The following block does delta decoding. Code is duplicated here because this is faster. */
	if (x >> 63 == 0) {
		int z = rld_ddec_table[x >> 55];
		if (z < 0) {
			z = ~z, w = z & 0xf, y = z >> 4;
		} else { // when z == 0, w = 0
			int a = z>>4, b = z&0xf;
			w = a + b, y = x << b >> (64 - a) | 1u << a;
		}
	} else w = y = 1;
	if (w == 0) return 0;
	y = y << e->abits | x << w >> (64 - e->abits);
	w += e->abits;
	if (e->r > w) e->r -= w;
	else ++e->p, e->r = 64 + e->r - w;
	return y;
}

static inline uint32_t rld_dec0(rld_t *e)
{
	int y = 0, w, l;
	uint64_t x;
	x = e->p[0] << (64 - e->r) | (e->p < e->stail && e->r < 64? e->p[1] >> e->r : 0);
	if (x>>63 == 0) {
		//w = 0x3333333355557790ll>>(x>>59<<2)&0xf;
		if (x>>62 == 0) {
			if (x>>61 == 0) {
				if ((w = 0x55555555777799b0ll>>(x>>58<<2)&0xf) == 0) return 0;
			} else w = 5;
		} else w = 3;
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

#endif
