#include <stdlib.h>
#include <assert.h>
#include "rle6.h"

#define NDEBUG

rle6_t *rle6_enc_init()
{
	rle6_t *r;
	r = (rle6_t*)calloc(1, sizeof(rle6_t));
	r->n = 1;
	r->z = (uint8_t**)malloc(sizeof(void*));
	r->lhead = r->shead = r->z[0] = (uint8_t*)calloc(RLE6_LSIZE, 1);
	r->p = r->shead + 48;
	r->stail = r->shead + RLE6_SSIZE;
	return r;
}

int rle6_enc(rle6_t *r, int l, int c)
{
	int i, w;
	assert(c < 6);
	assert(l < 1<<23);
	assert(r->p > r->shead && r->p < r->stail);
	w = l < 32? 1 : l < 1024? 2 : 4;
	if (r->p + w >= r->stail) {
		uint64_t *q = (uint64_t*)r->shead;
		for (i = 0; i < 6; ++i) q[i] = r->cnt[i];
		*r->p = 0xff;
		if (r->p + w - r->lhead >= RLE6_LSIZE) {
			++r->n;
			r->z = (uint8_t**)realloc(r->z, r->n * sizeof(void*));
			r->lhead = r->shead = r->z[r->n - 1] = (uint8_t*)calloc(RLE6_LSIZE, 1);
		} else r->shead += RLE6_SSIZE;
		r->p = r->shead + 48;
		r->stail = r->shead + RLE6_SSIZE;
	}
	if (w == 1) *r->p++ = c<<5 | l;
	else if (w == 2) {
		*r->p++ = 6<<5 | c<<2 | (l&0x3);
		*r->p++ = l>>2&0xff;
	} else {
		*r->p++ = 6<<5 | 6<<2 | (c&0x3);
		*r->p++ = (c>>2)<<7 | l&0x7f;
		*r->p++ = l>>7&0xff;
		*r->p++ = l>>15;
	}
	++r->cnt[c];
	return 0;
}

uint64_t rle6_enc_finish(rle6_t *r)
{
	int i;
	uint64_t *q = (uint64_t*)r->shead;
	for (i = 0; i < 6; ++i) q[i] = r->cnt[i];
	*r->p++ = 0xff;
	return (uint64_t)(r->n - 1) * RLE6_LSIZE + (r->p - r->lhead);
}

#include <stdio.h>
int main()
{
	int i, n = 20000;
	rle6_t *r = rle6_enc_init();
	for (i = 1; i < n; ++i)
		rle6_enc(r, i, 0);
	rle6_enc_finish(r);
	rle6_dec_init(r, 0);
	for (i = 1; i < n; ++i)
		printf("%d\t%d\n", i, rle6_dec(r)>>3);
	return 0;
}
