#include <stdlib.h>
#include <stdint.h>
#include "rld.h"

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2(uint32_t v)
{
	register uint32_t t, tt;
	if (tt = v >> 16) return (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v >> 8) ? 8 + LogTable256[t] : LogTable256[v];
}

static inline uint32_t delta_enc1(uint32_t x, int *width)
{
	int y = ilog2(x);
	int z = ilog2(y + 1);
	*width = (z<<1) + 1 + y;
	return x^(1<<y) | (y+1)<<y;
}

rldenc_t *rld_enc_init(int asize, int bbits)
{
	rldenc_t *e;
	e = calloc(1, sizeof(rldenc_t));
	e->z = malloc(sizeof(void*));
	e->z[0] = calloc(RLD_SUPBLK_SIZE, 8);
	e->p = e->z[0] + asize;
	e->cnt = calloc(asize, 8);
	e->b = ilog2(asize) + 1;
	e->r = 64;
	e->asize = asize;
	e->bbits = bbits;
	e->bmask = (1<<bbits) - 1;
	return r;
}

int rld_push(rldenc_t *e, int l, uint8_t c)
{
	int i, w;
	uint64_t x = delta_enc1(l, &w) << r->b | c;
	w += r->b;
	if (w > e->r) {
		if (((p - e->z[e->n] + 1) & e->bmask) == 0) { // jump to the next block
			int y = (p - e->z[e->n]) & e->bmask;
			for (i = 0; i < e->asize; ++i)
				e->z[e->n][y + i] = e->cnt[i];
			if (p - e->z[e->n] + 1 == RLD_SUPBLK_SIZE) { // allocate a new superblock
				++e->n;
				e->z = realloc(e->z, (e->n + 1) * sizeof(void*));
				p = e->z[e->n] = calloc(RLD_SUPBLK_SIZE, 8);
			} else ++p;
			e->p += e->asize;
			e->r = 64 - w;
			*e->p |= x << e->r;
		} else {
			w -= e->r;
			*e->p++ |= x >> w;
			*e->p = x << (e->r = 64 - w);
		}
	} else e->r -= w, *e->p |= x << e->r;
	++e->cnt[c];
	return 0;
}

uint64_t rld_finish(rldenc_t *e)
{
	int i, y = (p - e->z[e->n]) & e->bmask;
	for (i = 0; i < e->asize; ++i)
		e->z[e->n][y + i] = e->cnt[i];
	return (((uint64_t)(e->n + 1) * RLD_SUPBLK_SIZE) + (p - e->z[e->n])) * 64 + (64 - e->r);
}

/*
int main(int argc, char *argv[])
{
	int w, x;
	if (argc > 1) {
		x = delta_enc1(atoi(argv[1]), &w);
		printf("%x, %d\n", x, w);
	}
	return 0;
}
*/
