#ifndef BIT3_H
#define BIT3_H

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#define B3_SSIZE 32
#define B3_FSIZE 16
#define B3_SFSIZE (B3_SSIZE * B3_FSIZE)
#define B3_LSIZE (1<<25)

typedef struct {
	int n;
	int64_t n_bytes;
	uint16_t **z;
	int64_t n_frames;
	uint64_t *frame;
	int64_t *cnt, *mcnt;
	uint64_t *cnttab;
} bit3_t;

typedef struct {
	int r;
	uint16_t *p, **i;
} b3itr_t;

bit3_t *b3_init();
void b3_destroy(bit3_t *b3);
int64_t b3_enc_finish(bit3_t *b3, b3itr_t *itr);
int b3_dump(const bit3_t *b3, const char *fn);
bit3_t *b3_restore(FILE *fp);

static inline int b3_enc(bit3_t *b3, b3itr_t *itr, int64_t l, int c)
{
	int64_t i;
	for (i = 0; i < l; ++i) {
		*itr->p |= c<<itr->r;
		if (itr->r == 12) {
			itr->r = 0; ++itr->p;
			if (itr->p - *itr->i == B3_LSIZE) {
				++b3->n;
				b3->z = realloc(b3->z, sizeof(void*) * b3->n);
				b3->z[b3->n - 1] = calloc(B3_LSIZE, 2);
				itr->i = b3->z + b3->n - 1;
				itr->p = *itr->i;
			}
		} else itr->r += 3;
	}
	return 0;
}

static inline int b3_dec0(const bit3_t *b3, b3itr_t *itr)
{
	int c = *itr->p>>itr->r & 7;
	if (itr->r == 12) {
		itr->r = 0; ++itr->p;
		if (itr->p - *itr->i == B3_LSIZE) itr->p = *++itr->i;
	} else itr->r += 3;
	assert((itr->p - *itr->i) * 2 < b3->n_bytes);
	return c;
}

static inline int b3_dec(const bit3_t *b3, b3itr_t *itr, int *c, int is_free)
{
	*c = b3_dec0(b3, itr);
	return *c == 7? -1 : 1;
}

#endif
