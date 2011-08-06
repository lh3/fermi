#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#define RLD_LBITS 23
#define RLD_LSIZE (1<<RLD_LBITS)
#define RLD_LMASK (RLD_LSIZE - 1)

typedef struct {
	// static members (initialized in the constructor)
	int16_t asize, asize1; // alphabet size
	int abits; // bits required to store a symbol
	int sbits; // bits per small block
	int ssize; // ssize = 1<<sbits
	int8_t o0[2];
	// dynamic members
	int n; // number of blocks (unchanged in decoding)
	int r; // bits remaining in the last 64-bit integer
	uint64_t n_bytes; // total number of bits (unchanged in decoding)
	uint64_t **z; // the actual data (unchanged in decoding)
	uint64_t *cnt, *mcnt; // after enc_finish, cnt keeps the accumulative count and mcnt keeps the marginal
	uint64_t *p, *lhead, *shead, *stail;
} rld_t;

typedef struct {
	int b;
	uint64_t n, *s;
} rldidx_t;

#ifdef __cplusplus
extern "C" {
#endif

	rld_t *rld_enc_init(int asize, int bbits);
	int rld_enc(rld_t *e, int l, uint8_t c);
	uint64_t rld_enc_finish(rld_t *e);
	uint64_t rld_rawlen(const rld_t *e);
	rldidx_t *rld_index(const rld_t *e);

	int rld_dump(const rld_t *e, const char *fn);
	rld_t *rld_restore(const char *fn);
	void rld_destroy(rld_t *e);

#ifdef __cplusplus
}
#endif

#define rld_last_blk(e) ((e)->n_bytes>>3>>(e)->sbits<<(e)->sbits)
#define rld_seek_blk(e, k) ((e)->z[(k)>>RLD_LBITS] + ((k)&RLD_LMASK))

#define rld_size_bit(x) ((x)>>31&1)

static inline void rld_dec_init(rld_t *e, uint64_t k)
{
	e->lhead = e->z[k>>RLD_LBITS];
	e->shead = e->lhead + k%RLD_LSIZE;
	e->stail = e->shead + e->ssize - 1;
	e->p = e->shead + e->o0[rld_size_bit(*e->shead)];
	e->r = 64;
}

static inline int rld_dec0(rld_t *e, int *c)
{
	int y = 0, w, l;
	uint64_t x;
	assert(e->p <= e->stail);
	x = e->p[0] << (64 - e->r) | (e->p < e->stail && e->r < 64? e->p[1] >> e->r : 0);
	if (x>>63 == 0) {
		if ((w = 0x333333335555779bll>>(x>>59<<2)&0xf) == 0xb && x>>58 == 0) return 0;
		l = (x >> (64 - w)) - 1;
		y = x << w >> (64 - l) | 1u << l;
		w += l;
	} else w = y = 1;
	*c = x << w >> (64 - e->abits);
	w += e->abits;
	if (e->r > w) e->r -= w;
	else ++e->p, e->r = 64 + e->r - w;
	return y;
}

static inline int rld_dec(rld_t *e, int *_c)
{
	int l = rld_dec0(e, _c);
	if (l == 0) {
		uint64_t last = rld_last_blk(e);
		if (e->p - e->lhead > RLD_LSIZE - e->ssize) e->shead = ++e->lhead;
		else e->shead += e->ssize;
		if (e->shead == rld_seek_blk(e, last)) return -1;
		e->p = e->shead + e->o0[rld_size_bit(*e->shead)];
		e->stail = e->shead + e->ssize - 1;
		e->r = 64;
		return rld_dec0(e, _c);
	} else return l;
}

static inline uint64_t *rld_locate_blk(rld_t *e, const rldidx_t *r, uint64_t k, uint64_t *cnt, uint64_t *sum)
{
	int j;
	uint64_t *q, *z = r->s + (k>>r->b) * e->asize1;
	e->lhead = e->z[*z>>RLD_LBITS];
	q = e->p = e->lhead + (*z&RLD_LMASK);
	for (j = 1, *sum = 0; j < e->asize1; ++j) *sum += (cnt[j-1] = z[j]);
	while (1) { // seek to the small block
		uint64_t c = 0;
		if (q - e->lhead == RLD_LSIZE) q = ++e->lhead;
		else q += e->ssize;
		c = rld_size_bit(*q)? *((uint32_t*)q)&0x7fffffff : *(uint16_t*)q;
		if (*sum + c > k) break;
		if (rld_size_bit(*q)) {
			uint32_t *p = (uint32_t*)q;
			for (j = 0; j < e->asize; ++j) cnt[j] += p[j+1];
		} else {
			uint16_t *p = (uint16_t*)q;
			for (j = 0; j < e->asize; ++j) cnt[j] += p[j+1];
		}
		*sum += c;
		e->p = q;
	}
	e->shead = e->p;
	e->stail = e->shead + e->ssize - 1;
	e->p += e->o0[*e->shead>>63];
	e->r = 64;
	return q;
}

static inline uint64_t rld_rank11(rld_t *e, const rldidx_t *r, uint64_t k, int c)
{
	uint64_t y, z, *cnt;
	cnt = alloca(e->asize1 * 8);
	rld_locate_blk(e, r, k, cnt, &z);
	y = cnt[c];
	++k; // because k is the coordinate but not length
	while (1) {
		int a = -1, l;
		l = rld_dec0(e, &a);
		if (z + l >= k) return y + (a == c? k - z: 0);
		z += l;
		if (a == c) y += l;
	}
}

static inline void rld_rank1a(rld_t *e, const rldidx_t *r, uint64_t k, uint64_t *ok)
{
	uint64_t z;
	rld_locate_blk(e, r, k, ok, &z);
	++k; // because k is the coordinate but not length
	while (1) {
		int a = -1, l;
		l = rld_dec0(e, &a);
		if (z + l >= k) {
			ok[a] += k - z;
			break;
		}
		z += l; ok[a] += l;
	}
}

static inline void rld_rank21(rld_t *e, const rldidx_t *r, uint64_t k, uint64_t l, int c, uint64_t *ok, uint64_t *ol) // FIXME: can be faster
{
	*ok = rld_rank11(e, r, k, c);
	*ol = rld_rank11(e, r, l, c);
}

static inline void rld_rank2a(rld_t *e, const rldidx_t *r, uint64_t k, uint64_t l, uint64_t *ok, uint64_t *ol) // FIXME: can be faster
{
	rld_rank1a(e, r, k, ok);
	rld_rank1a(e, r, l, ol);
}

#endif
