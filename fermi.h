#ifndef FERMI_H
#define FERMI_H

#include <stdint.h>
#include <stdlib.h>

#define FERMI_VERSION "1.1-r751-beta"

#define FM_MASK30 0x3fffffff

extern int fm_verbose;

typedef struct {
	uint64_t x[3]; // 0: start of the interval, backward; 1: forward; 2: size of the interval
	uint64_t info;
} fmintv_t;

typedef struct { size_t n, m; int32_t *a; } fm32s_v;
typedef struct { size_t n, m; uint32_t *a; } fm32_v;
typedef struct { size_t n, m; uint64_t *a; } fm64_v;
typedef struct { size_t n, m; fmintv_t *a; } fmintv_v;

struct __rld_t; // defined in rld.h
struct __mog_t; // defined in mog.h

typedef struct {
	int w, min_occ, keep_bad, is_paired, trim_l, step;
	float max_corr;
} fmecopt_t;

typedef struct {
	int pr_links, min_supp;
	double avg, std, a_thres, p_thres;
} fmscafopt_t;

#ifndef KINT_DEF
#define KINT_DEF
typedef struct { uint64_t x, y; } ku128_t;
typedef struct { size_t n, m; uint64_t *a; } ku64_v;
typedef struct { size_t n, m; ku128_t *a; } ku128_v;
#endif

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

extern unsigned char seq_nt6_table[128];

#ifdef __cplusplus
extern "C" {
#endif

	int fm_bwtgen(int asize, int64_t l, uint8_t *s);
	struct __rld_t *fm_bwtenc(int asize, int sbits, int64_t l, const uint8_t *s);
	struct __rld_t *fm_append(struct __rld_t *e0, int len, const uint8_t *T);
	struct __rld_t *fm_build(struct __rld_t *e0, int asize, int sbits, int64_t l, uint8_t *s);
	struct __rld_t *fm6_build2(int64_t l, const char *s);

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
	 * @param s  output string (reversed)
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

	fmintv_t *fm6_traverse(const struct __rld_t *e, int depth);

	int fm6_smem1(const struct __rld_t *e, int len, const uint8_t *q, int x, fmintv_v *mem, int self_match);
	int fm6_smem(const struct __rld_t *e, int len, const uint8_t *q, fmintv_v *mem, int self_match);
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

	struct __rld_t *fm_sub(struct __rld_t *e, const uint64_t *sub, int n_threads, int is_comp);

	int64_t fm6_api_readseq(const char *fn, char **_seq, char **_qual);
	void fm6_api_writeseq(int64_t l, char *seq, char *qual);
	int fm6_api_seqlen(int64_t l, const char *seq, double quantile);
	int fm6_api_correct(int kmer, int64_t l, char *_seq, char *_qual);
	struct __mog_t *fm6_api_unitig(int min_match, int64_t l, char *seq);

#ifdef __cplusplus
}
#endif

#endif
