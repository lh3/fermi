#ifndef FM_PRIV_H
#define FM_PRIV_H

#include "rld.h"
#include "fermi.h"
#include "mag.h"
#include "utils.h"

// SA construction

int ksa_sa(const unsigned char *T, int *SA, int n, int k);
int ksa_bwt(unsigned char *T, int n, int k);
int ksa_bwt64(unsigned char *T, int64_t n, int k);

// Basic algorithms: sorting and heap

void ks_introsort_uint64_t(size_t n, uint64_t *a);
void ks_introsort_128x(size_t n, ku128_t *a);
void ks_introsort_128y(size_t n, ku128_t *a);
void ks_heapup_uint64_t(size_t n, uint64_t *a);
void ks_heapdown_uint64_t(size_t i, size_t n, uint64_t *a);
void ks_heapmake_uint64_t(size_t n, uint64_t *a);
void ks_heapup_128y(size_t n, ku128_t *a);
void ks_heapdown_128y(size_t i, size_t n, ku128_t *a);
void ks_heapmake_128y(size_t n, ku128_t *a);

// Basic sequence manipulation

void seq_char2nt6(int l, unsigned char *s);
void seq_reverse(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

// Portals to CLI algorithms

uint64_t *fm6_seqsort(const rld_t *e, int n_threads);
int fm6_unitig(const struct __rld_t *e, int min_match, int n_threads);
int fm6_ec_correct(const struct __rld_t *e, fmecopt_t *opt, const char *fn, int n_threads, const char *fn_hash);
int fm6_remap(const char *fn, const rld_t *e, uint64_t *sorted, int skip, int min_pcv, int max_dist, int n_threads);
void mag_scaf_core(const rld_t *e, const char *fn, const fmscafopt_t *opt, int n_threads);

void fm_reverse_fmivec(fmintv_v *p);

uint64_t fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s, fmintv_t *k2, int *contained);

void fm6_unpack_rlo(const rld_t *e);

// For the second correction algorithm

typedef struct {
	int min_l, min_occ, min_occ_patch;
	int max_pen, max_d, len_factor, diff_factor;
	int qual_plus;
	double ratio_factor;
} fmec2opt_t;

typedef struct {
	uint8_t oq, cb[2], cq[2], cf;
} fmec2seq_t;

typedef struct {
	int max_len, max_matrix;
	uint8_t *seq;
	fmec2seq_t *s;
	int *matrix;
	fmintv_v tmp[2];
	fmsmem_v mem1, mem;
	fm32s_v f, b; // for dynamic programming
} fmec2aux_t;

void fmc_opt_init(fmec2opt_t *opt);
void fmc_aux_destroy(fmec2aux_t *a);
void fmc_ec_core(const fmec2opt_t *opt, const rld_t *e, fmec2aux_t *aux, int l_seq, char *seq, char *qual);

#endif
