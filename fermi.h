#ifndef FERMI_H
#define FERMI_H

#include <stdint.h>
#include <stdlib.h>

#define FERMI_VERSION "0.0-dev (r341)"

#define FM_MASK30 0x3fffffff

extern int fm_verbose;

typedef struct {
	uint64_t x[3]; // 0: start of the interval, backward; 1: forward; 2: size of the interval
	uint64_t info;
} fmintv_t;

typedef struct {
	uint64_t x, y;
} fm128_t;

typedef struct { size_t n, m; fmintv_t *a; } fmintv_v;
typedef struct { size_t n, m; uint64_t *a; } fm64_v;
typedef struct { size_t n, m; fm128_t  *a; } fm128_v;
struct __rld_t; // defined in rld.h

typedef struct {
	int T, t, w, max_pre_mm;
	float cov, err;
} fmecopt_t;

typedef struct {
	float min_tip_cov, min_br_ratio;
	int min_tip_len, min_br_width;
	int check;
} fmclnopt_t;

typedef struct {
	uint64_t k[2];
	fm128_v nei[2];
	int l;
	float avg_cov;
	char *seq, *cov;
} fmnode_t;

typedef struct { size_t n, m; fmnode_t *a; } fmnode_v;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t { // implemented in kstring.h
	uint32_t l, m;
	char *s;
} kstring_t;
#endif

// complement of a nucleotide
#define fm6_comp(a) ((a) >= 1 && (a) <= 4? 5 - (a) : (a))
#define fm6_set_intv(e, c, ik) ((ik).x[0] = (e)->cnt[(int)(c)], (ik).x[2] = (e)->cnt[(int)(c)+1] - (e)->cnt[(int)(c)], (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

#ifdef __cplusplus
extern "C" {
#endif

	int fm_bwtgen(int asize, int64_t l, uint8_t *s);
	struct __rld_t *fm_bwtenc(int asize, int sbits, int64_t l, const uint8_t *s);
	struct __rld_t *fm_append(struct __rld_t *e0, int len, const uint8_t *T);
	struct __rld_t *fm_build(struct __rld_t *e0, int asize, int sbits, int64_t l, uint8_t *s);

	/**
	 * Backward search for a generic FM-Index
	 *
	 * @param e       FM-Index
	 * @param len     length of the input string
	 * @param str     input string
	 * @param sa_beg  the start of the final SA interval
	 * @param sa_end  the end of the interval
	 * 
	 * @return        equal to (*sa_end - *sa_end - 1); zero if no match
	 */
	uint64_t fm_backward_search(const struct __rld_t *e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end);

	/**
	 * Retrieve the x-th string from a generic FM-Index
	 *
	 * @param e  FM-Index
	 * @param x  string to retrieve (x >= 0)
	 * @param s  output string
	 */
	int64_t fm_retrieve(const struct __rld_t *e, uint64_t x, kstring_t *s);

	/**
	 * Extend a string in either forward or backward direction
	 *
	 * @param e        DNA FM-Index
	 * @param ik       input SA interval
	 * @param ok       output SA intervals; one for each symbol between 0 and 5
	 * @param is_back  true is backward (right-to-left); otherwise forward (left-to-right)
	 */
	int fm6_extend(const struct __rld_t *e, const fmintv_t *ik, fmintv_t ok[6], int is_back);
	int fm6_extend0(const struct __rld_t *e, const fmintv_t *ik, fmintv_t *ok0, int is_back);

	int fm6_smem1(const struct __rld_t *e, int len, const uint8_t *q, int x, fmintv_v *mem);
	int fm6_smem(const struct __rld_t *e, int len, const uint8_t *q, fmintv_v *mem);
	int fm6_write_smem(const struct __rld_t *e, const fmintv_t *a, kstring_t *s);

	/**
	 * Merge two generic FM-Indexes
	 *
	 * @param e0   first FM-Index
	 * @param e1   second FM-Index
	 * @param n_threads  number of threads to use
	 *
	 * @return     output FM-Index
	 */
	struct __rld_t *fm_merge(struct __rld_t *e0, struct __rld_t *e1, int n_threads);

	int fm6_unitig(const struct __rld_t *e, int min_match, int n_threads);

	void fm_ec_genpar(int64_t n, int l, double cov, double p, int *_w, int *_T);
	int fm6_ec_correct(const struct __rld_t *e, const fmecopt_t *opt, int n_threads);

	fmnode_v *msg_read(const char *fn);
	void msg_print(const fmnode_v *nodes);
	void msg_clean(fmnode_v *nodes, const fmclnopt_t *opt);

#ifdef __cplusplus
}
#endif

#endif
