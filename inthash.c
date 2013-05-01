#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include "inthash.h"

#ifndef kcalloc
#define kcalloc(type, p, cnt, size) if (((p) = (type)calloc((cnt), (size))) == 0) fprintf(stderr, "[E::%s] fail to allocate %ld bytes\n", __func__, (long)(size))
#endif

// each uint32_t integer is separated to: "rest_key << (b_idx+b_val) | idx << b_idx | val << b_val"
static inline const uint32_t *ih_get_core(int b_tbl, const uint32_t *h, int b_idx, int b_val, uint64_t x)
{ // note that we are using a table size of power of 2 and triangular steps; we always traverse each bucket (see wiki)
	uint64_t mask = (1ULL<<b_tbl) - 1, z = x&mask, z0 = z;
	uint32_t y = x >> b_tbl;
	int i = 0;
	while (h[z] && h[z]>>b_val != (y<<b_idx|(i+1))) {
		z = (z + (++i)) & mask;
		if (z == z0 || i+1 >= 1U<<b_idx) return 0;
	}
	return &h[z];
}

static inline uint32_t *ih_put_core(int b_tbl, uint32_t *h, int b_idx, int b_val, uint64_t x, int *present)
{
	uint64_t mask = (1ULL<<b_tbl) - 1, z = x&mask, z0 = z;
	uint32_t y = x >> b_tbl;
	int i = 0;
	*present = 0;
	while (h[z] && h[z]>>b_val != (y<<b_idx|(i+1))) {
		z = (z + (++i)) & mask;
		if (z == z0 || i+1 >= 1U<<b_idx) return 0;
	}
	if (h[z] == 0) h[z] = y<<(b_val+b_idx) | (i+1)<<b_val;
	else *present = 1;
	return &h[z];
}

static uint32_t *ih_dbl_core(int b_tbl, uint32_t *h0, int b_idx, int b_val)
{
	uint64_t mask = (1ULL<<b_tbl) - 1, z, old_size = 1ULL<<b_tbl;
	uint32_t *h, lmask = (1U<<b_idx) - 1, vmask = (1U<<b_val) - 1;
	kcalloc(uint32_t*, h, old_size<<1, 4);
	for (z = 0; z < old_size; ++z) {
		uint64_t x, i;
		uint32_t *p;
		int present;
		if (h0[z] == 0) continue;
		i = h0[z]>>b_val & lmask;
		x = (z - (i*(i-1)>>1)) & mask;
		x = h0[z] >> (b_idx+b_val) << b_tbl | x;
		p = ih_put_core(b_tbl+1, h, b_idx, b_val, x, &present);
		*p |= h0[z] & vmask;
	}
	free(h0);
	return h;
}

inthash_t *ih_init(int b_key, int b_idx, int b_val)
{
	inthash_t *ih;
	int b_tbl = b_key + b_idx + b_val - 32;
	if (b_tbl >= 32) return 0;
	if (b_tbl < 4) b_tbl = 4;
	ih = calloc(1, sizeof(inthash_t));
	ih->b_val = b_val; ih->b_key = b_key; ih->b_tbl = b_tbl; ih->b_idx = b_idx;
	ih->h = calloc(1ULL<<b_tbl, 4);
	return ih;
}

uint32_t *ih_put(inthash_t *ih, uint64_t x)
{
	uint32_t *p;
	int b_tbl, present;
	while (__sync_lock_test_and_set(&ih->lock, 1)) while (ih->lock); // spin lock
	b_tbl = ih->b_tbl;
	if (ih->n_elem > (1ULL<<ih->b_tbl>>1) + (1ULL<<ih->b_tbl>>2)) // rehashing if 3/4 full
		ih->h = ih_dbl_core(b_tbl++, ih->h, ih->b_idx, ih->b_val);
	if ((p = ih_put_core(b_tbl, ih->h, ih->b_idx, ih->b_val, x, &present)) == 0) {
		ih->h = ih_dbl_core(b_tbl++, ih->h, ih->b_idx, ih->b_val);
		p = ih_put_core(b_tbl, ih->h, ih->b_idx, ih->b_val, x, &present);
	}
	ih->n_elem += !present;
	ih->b_tbl = b_tbl;
	__sync_lock_release(&ih->lock);
	return p;
}

const uint32_t *ih_get(const inthash_t *ih, uint64_t x)
{
	return ih_get_core(ih->b_tbl, ih->h, ih->b_idx, ih->b_val, x);
}
/*
int main()
{
	int i, N = 10000, mask = 0x1ff;
	inthash_t *ih;
	srand48(11);
	ih = ih_init(20, 6, 8);
	fprintf(stderr, "Initial table bits: %d\n", ih->b_tbl);
	for (i = 0; i < N; ++i) {
		uint64_t x = lrand48() & mask;
		ih_put(ih, x);
		printf("I\t%lld\n", x);
	}
	for (i = 0; i < 1<<ih->b_tbl; ++i)
		if (ih->h[i]) printf("O\t%llu\n", ih_dec(ih, &ih->h[i]));
	fprintf(stderr, "Final table bits: %d (%u elements)\n", ih->b_tbl, ih->n_elem);
	ih_destroy(ih);
	return 0;
}
*/
