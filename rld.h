#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>
#define RLD_SUPBLK_BIT  23
#define RLD_SUPBLK_SIZE (1<<RLD_SUPBLK_BIT)

typedef struct {
	int asize; // alphabet size
	int bbits; // bits per small block
	int bmask;
	int n; // current block
	int b; // bits required to store a symbol
	int r; // bits remaining in the last 64-bit integer
	uint64_t *cnt;
	uint64_t **z;
	uint64_t *p;
} rldenc_t;

rldenc_t *rld_enc_init(int asize, int bbits);
int rld_push(rldenc_t *e, int l, uint8_t c);
uint64_t rld_finish(rldenc_t *e);

#endif
