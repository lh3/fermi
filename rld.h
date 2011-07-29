#ifndef RLDELTA_H
#define RLDELTA_H

#include <stdint.h>
#define RLD_SUPBLK_BIT  23
#define RLD_SUPBLK_SIZE (1<<RLD_SUPBLK_BIT)
#define rld_blk_size(bbits, balpha) (((1<<(bbits)) << 6) / ((balpha) + 1))

typedef struct {
	// static members
	int asize; // alphabet size
	int bbits; // bits per small block
	int bmask;
	int b; // bits required to store a symbol
	// dynamic members in encoding
	int n; // number of blocks
	// dynamic members in decoding
	int r; // bits remaining in the last 64-bit integer
	uint64_t *cnt;
	uint64_t **z;
	uint64_t *p, *head;
} rld_t;

#ifdef __cplusplus
extern "C" {
#endif

	rld_t *rld_enc_init(int asize, int bbits);
	int rld_push(rld_t *e, int l, uint8_t c);
	uint64_t rld_finish(rld_t *e);

#ifdef __cplusplus
}
#endif

#endif
