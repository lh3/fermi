#ifndef FM_MOG_H
#define FM_MOG_H

#include <stdint.h>
#include <stdlib.h>
#include "kstring.h"

#define MOG_F_DROP_TIP0  0x1
#define MOG_F_READ_TAG   0x2
#define MOG_F_READnMERGE 0x4
#define MOG_F_CLEAN      0x10
#define MOG_F_AGGRESSIVE 0x20

typedef struct {
	int flag, max_arc, n_iter, min_ovlp, min_elen, min_ensr, min_insr;
	float min_dratio0, min_dratio1, a_thres;
} mogopt_t;

#ifndef KINT_DEF
#define KINT_DEF
typedef struct { uint64_t x, y; } ku128_t;
typedef struct { size_t n, m; uint64_t *a; } ku64_v;
typedef struct { size_t n, m; ku128_t *a; } ku128_v;
#endif

typedef struct {
	int len, nsr;    // length; number supporting reads
	uint32_t max_len;// allocated seq/cov size
	int aux[3];      // auxiliary information
	uint64_t k[2];   // bi-interval
	ku128_v nei[2];  // neighbors
	char *seq, *cov; // sequence and coverage
	void **ptr;      // additional information
} mogv_t;

typedef struct { size_t n, m; mogv_t *a; } mogv_v;

typedef struct {
	mogv_v v;
	float rdist;  // read distance
	int min_ovlp; // minimum overlap seen from the graph
	void *h;
} mog_t;

#ifdef __cplusplus
extern "C" {
#endif
	mogopt_t *mog_init_opt(void);
	void mog_g_clean(mog_t *g, const mogopt_t *opt);

	void mog_g_destroy(mog_t *g);
	mog_t *mog_g_read(const char *fn, const mogopt_t *opt);
	void mog_v_write(const mogv_t *p, kstring_t *out);
	void mog_g_print(const mog_t *g);
	void mog_g_merge(mog_t *g);
	double mog_cal_rdist(mog_t *g);

	void mog_v_copy_to_empty(mogv_t *dst, const mogv_t *src); // NB: memory leak if dst is allocated

#ifdef __cplusplus
}
#endif

#endif
