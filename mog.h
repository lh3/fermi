#ifndef FM_MOG_H
#define FM_MOG_H

#include <stdint.h>
#include <stdlib.h>
#include "kstring.h"

#define MOG_F_DROP_TIP 0x1
#define MOG_F_READ_TAG 0x2

typedef struct {
	uint64_t x, y;
} mog128_t;

typedef struct {
	int flag, max_arc, min_el;
	float diff_ratio;
} mogopt_t;

typedef struct { size_t n, m; mog128_t *a; } mog128_v;

typedef struct {
	int len, nsr;    // length; number supporting reads
	int max_len;
	uint64_t k[2];   // bi-interval
	mog128_v nei[2]; // neighbors
	char *seq, *cov; // sequence and coverage
	int aux[2];      // auxiliary information
	void **ptr;      // additional information
} mognode_t;

typedef struct { size_t n, m; mognode_t *a; } mognode_v;

typedef struct {
	mognode_v v;
	void *h;
	int min_ovlp;
} mog_t;

void mog_write1(const mognode_t *p, kstring_t *out);

#endif
