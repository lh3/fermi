#ifndef RLE6_H
#define RLE6_H

#define RLE6_LSHIFT 24
#define RLE6_LSIZE  (1<<RLE6_LSHIFT)

#define RLE6_SSHIFT 8
#define RLE6_SSIZE  (1<<RLE6_SSHIFT)

#include <stdint.h>

typedef struct {
	int n;
	uint8_t **z;
	uint8_t *p, *lhead, *shead, *stail;
	uint64_t cnt[6];
} rle6_t;

#ifdef __cplusplus
extern "C" {
#endif

	rle6_t *rle6_enc_init();
	int rle6_enc(rle6_t *r, int l, int c);
	uint64_t rle6_enc_finish(rle6_t *r);

#ifdef __cplusplus
}
#endif

static inline int rle6_dec_init(rle6_t *r, uint64_t x)
{
	x = x >> RLE6_SSHIFT << RLE6_SSHIFT;
	r->lhead = r->z[x>>RLE6_LSHIFT];
	r->shead = r->lhead + x%RLE6_LSIZE;
	r->stail = r->shead + RLE6_SSIZE;
	r->p = r->shead + 48;
	return 0;
}

static inline uint32_t rle6_dec0(rle6_t *r)
{
	uint32_t c = *r->p >> 5;
	if (c < 6) {
		c |= (*r->p&0x1f) << 3;
		++r->p;
	} else if (c == 6) {
		c = *r->p>>2 & 0x7;
		if (c != 6) {
			c |= ((uint32_t)r->p[0]&0x3)<<3 | (uint32_t)r->p[1]<<5;
			r->p += 2;
		} else {
			c = (r->p[0]&0x3) | r->p[1]>>7<<2 | ((uint32_t)r->p[1]&0x7f) << 3 | (uint32_t)r->p[2] << 10 | (uint32_t)r->p[3] << 18;
			r->p += 4;
		}
	} else return 0;
	return c;
}

static inline uint32_t rle6_dec(rle6_t *r)
{
	uint32_t c = rle6_dec0(r);
	if (c == 0) {
		if (r->p - r->lhead > RLE6_LSIZE - RLE6_SSIZE) {
			++r->lhead;
			r->shead = r->lhead;
		} else r->shead += RLE6_SSIZE;
		r->p = r->shead + 48;
		return rle6_dec0(r);
	} else return c;
}

#endif
