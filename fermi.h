#ifndef FERMI_H
#define FERMI_H

#include <stdint.h>

#define FERMI_VERSION "0.0-dev (r212)"

extern int fm_verbose;

typedef struct {
	uint64_t x[3]; // 0: start of the interval, backward; 1: forward; 2: size of the interval
} fmintv_t;

struct __rld_t; // defined in rld.h

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t { // implemented in kstring.h
	uint32_t l, m;
	char *s;
} kstring_t;
#endif

// complement of a nucleotide
#define fm6_comp(a) ((a) >= 1 && (a) <= 4? 5 - (a) : (a))

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
	void fm_retrieve(const struct __rld_t *e, uint64_t x, kstring_t *s);

	/** Similar to {@link #fm_retrieve()} but working for DNA FM-Index only */
	void fm6_retrieve(const struct __rld_t *e, uint64_t x, kstring_t *s);

	/**
	 * Extend a string in either forward or backward direction
	 *
	 * @param e        DNA FM-Index
	 * @param ik       input SA interval
	 * @param ok       output SA intervals; one for each symbol between 0 and 5
	 * @param is_back  true is backward (right-to-left); otherwise forward (left-to-right)
	 */
	int fm6_extend(const struct __rld_t *e, const fmintv_t *ik, fmintv_t ok[6], int is_back);

	/**
	 * Find the exact match in the left-to-right direction.
	 *
	 * Matching continues if sequences in the index have minimum {min} exact matches.
	 *
	 * @param e    DNA FM-Index
	 * @param min  minimum end-to-end overlaps between sequences in the index
	 * @param len  length of the input string (DNA sequence)
	 * @param seq  input string
	 *
	 * @return     length of the maximal match
	 */
	int fm6_search_overlap(const struct __rld_t *e, int min, int len, const uint8_t *seq, int is_back);

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

#ifdef __cplusplus
}
#endif

#endif
