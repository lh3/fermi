#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#define RLD_LBITS 23
#define RLD_LSIZE (1<<RLD_LBITS)
#define RLD_LMASK (RLD_LSIZE - 1)

typedef struct {
	int r; // bits remained in the last 64-bit integer
	uint64_t *p, *shead, *stail, **i;
} rlditr_t;

typedef struct {
	// initialized in the constructor
	uint8_t asize, asize1; // alphabet size; asize1=asize+1
	int8_t abits; // bits required to store a symbol
	int8_t sbits; // bits per small block
	int8_t ibits; // modified during indexing; here for a better alignment
	int8_t offset0[2];
	int ssize; // ssize = 1<<sbits
	// modified during encoding
	int n; // number of blocks (unchanged in decoding)
	uint64_t n_bytes; // total number of bits (unchanged in decoding)
	uint64_t **z; // the actual data (unchanged in decoding)
	uint64_t *cnt, *mcnt; // after enc_finish, cnt keeps the accumulative count and mcnt keeps the marginal
	// modified during indexing
	uint64_t n_frames;
	uint64_t *frame;
} rld_t;

#ifdef __cplusplus
extern "C" {
#endif

	rld_t *rld_init(int asize, int bbits);
	int rld_enc(rld_t *e, rlditr_t *itr, int l, uint8_t c);
	uint64_t rld_enc_finish(rld_t *e, rlditr_t *itr);

	int rld_dump(const rld_t *e, const char *fn);
	rld_t *rld_restore(const char *fn);
	void rld_destroy(rld_t *e);

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
	itr->r = 64;
}

static inline int rld_dec0(const rld_t *e, rlditr_t *itr, int *c)
{
	int y = 0, w, l;
	uint64_t x;
	assert(itr->p <= itr->stail);
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

static inline int rld_dec(const rld_t *e, rlditr_t *itr, int *_c)
{
	int l = rld_dec0(e, itr, _c);
	if (l == 0) {
		uint64_t last = rld_last_blk(e);
		if (itr->p - *itr->i > RLD_LSIZE - e->ssize) itr->shead = *++itr->i;
		else itr->shead += e->ssize;
		if (itr->shead == rld_seek_blk(e, last)) return -1;
		itr->p = itr->shead + e->offset0[rld_size_bit(*itr->shead)];
		itr->stail = itr->shead + e->ssize - 1;
		itr->r = 64;
		return rld_dec0(e, itr, _c);
	} else return l;
}

static inline uint64_t *rld_locate_blk(const rld_t *e, rlditr_t *itr, uint64_t k, uint64_t *cnt, uint64_t *sum)
{
	int j;
	uint64_t *q, *z = e->frame + (k>>e->ibits) * e->asize1;
	itr->i = e->z + (*z>>RLD_LBITS);
	q = itr->p = *itr->i + (*z&RLD_LMASK);
	for (j = 1, *sum = 0; j < e->asize1; ++j) *sum += (cnt[j-1] = z[j]);
	while (1) { // seek to the small block
		uint64_t c = 0;
		if (q - *itr->i == RLD_LSIZE) q = *++itr->i;
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
		itr->p = q;
	}
	itr->shead = itr->p;
	itr->stail = itr->shead + e->ssize - 1;
	itr->p += e->offset0[*itr->shead>>63];
	itr->r = 64;
	return q;
}

static inline uint64_t rld_rank11(const rld_t *e, uint64_t k, int c)
{
	uint64_t y, z, *cnt;
	rlditr_t itr;
	cnt = alloca(e->asize1 * 8);
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

static inline void rld_rank1a(const rld_t *e, uint64_t k, uint64_t *ok)
{
	uint64_t z;
	rlditr_t itr;
	rld_locate_blk(e, &itr, k, ok, &z);
	++k; // because k is the coordinate but not length
	while (1) {
		int a = -1, l;
		l = rld_dec0(e, &itr, &a);
		if (z + l >= k) {
			ok[a] += k - z;
			break;
		}
		z += l; ok[a] += l;
	}
}

static inline void rld_rank21(const rld_t *e, uint64_t k, uint64_t l, int c, uint64_t *ok, uint64_t *ol) // FIXME: can be faster
{
	*ok = rld_rank11(e, k, c);
	*ol = rld_rank11(e, l, c);
}

static inline void rld_rank2a(const rld_t *e, uint64_t k, uint64_t l, uint64_t *ok, uint64_t *ol) // FIXME: can be faster
{
	rld_rank1a(e, k, ok);
	rld_rank1a(e, l, ol);
}

#endif
