#ifndef FM_MOG_H
#define FM_MOG_H

#include <stdint.h>
#include <stdlib.h>
#include "kstring.h"

#define MOG_F_DROP_TIP0 0x1
#define MOG_F_READ_TAG  0x2
#define MOG_F_CLEAN     0x10

typedef struct {
	int flag, max_arc, min_el;
	float min_dratio0;
} mogopt_t;

#ifndef KINT_DEF
#define KINT_DEF
typedef struct { uint64_t x, y; } ku128_t;
typedef struct { size_t n, m; uint64_t *a; } ku64_v;
typedef struct { size_t n, m; ku128_t *a; } ku128_v;
#endif

typedef struct {
	int len, nsr;    // length; number supporting reads
	int max_len;
	uint64_t k[2];   // bi-interval
	ku128_v nei[2];  // neighbors
	char *seq, *cov; // sequence and coverage
	int aux[2];      // auxiliary information
	void **ptr;      // additional information
} mogv_t;

typedef struct { size_t n, m; mogv_t *a; } mogv_v;

typedef struct {
	mogv_v v;
	void *h;
	int min_ovlp;
} mog_t;

mogopt_t *mog_init_opt(void);
mog_t *mog_read(const char *fn, const mogopt_t *opt);
void mog_write1(const mogv_t *p, kstring_t *out);
void mog_print(const mogv_v *v);

#endif
