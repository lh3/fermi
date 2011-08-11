#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define RLD_LBITS  26
#define RLD_LSIZE  (1<<RLD_LBITS)
#define RLD_LMASK  (RLD_LSIZE - 1)

#define RLD_ABITS  3
#define RLD_ASIZE  6
#define RLD_ASIZE1 7
#define RLD_OFFSET 8

typedef struct {
	uint8_t *q, *shead, *stail, **i;
} rlditr_t;

typedef struct __rld_t {
	// initialized in the constructor
	int8_t sbits; // bits per small block
	int8_t ibits; // modified during indexing; here for a better alignment
	uint16_t ssize; // ssize = 1<<sbits
	// modified during encoding
	int n; // number of blocks
	uint64_t n_bytes; // total number of bits
	uint64_t *cnt, *mcnt; // after enc_finish, cnt keeps the accumulative count and mcnt keeps the marginal
	uint8_t **z; // the actual data
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
	int rld_enc(rld_t *e, rlditr_t *itr, int l, uint8_t c);
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

#define rld_last_blk(e) ((e)->n_bytes>>(e)->sbits<<(e)->sbits)
#define rld_seek_blk(e, k) ((e)->z[(k)>>RLD_LBITS] + ((k)&RLD_LMASK))

static inline void rld_itr_init(const rld_t *e, rlditr_t *itr, uint64_t k)
{
	itr->i = e->z + (k >> RLD_LBITS);
	itr->shead = *itr->i + k%RLD_LSIZE;
	itr->stail = itr->shead + e->ssize - 1;
	itr->q = itr->shead + RLD_OFFSET;
}

static inline int rld_dec0(const rld_t *r, rlditr_t *itr, int *c)
{
	*c = *itr->q >> 5;
	return *c == 7? 0 : *itr->q++ & 0x1f;
}

static inline int rld_dec(const rld_t *e, rlditr_t *itr, int *_c)
{
	int l = rld_dec0(e, itr, _c);
	if (l == 0) {
		uint64_t last = rld_last_blk(e);
		if (itr->q - *itr->i > RLD_LSIZE - e->ssize) itr->shead = *++itr->i;
		else itr->shead += e->ssize;
		if (itr->shead == rld_seek_blk(e, last)) return -1;
		itr->q = itr->shead + RLD_OFFSET;
		itr->stail = itr->shead + e->ssize - 1;
		return rld_dec0(e, itr, _c);
	} else return l;
}

#endif
