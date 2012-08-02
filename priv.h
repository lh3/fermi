#ifndef FM_PRIV_H
#define FM_PRIV_H

#include "rld.h"
#include "fermi.h"
#include "mag.h"
#include "utils.h"

int ksa_sa(const unsigned char *T, int *SA, int n, int k);
int ksa_bwt(unsigned char *T, int n, int k);
int ksa_bwt64(unsigned char *T, int64_t n, int k);

void ks_introsort_uint64_t(size_t n, uint64_t *a);
void ks_introsort_128x(size_t n, ku128_t *a);
void ks_introsort_128y(size_t n, ku128_t *a);
void ks_heapup_uint64_t(size_t n, uint64_t *a);
void ks_heapdown_uint64_t(size_t i, size_t n, uint64_t *a);
void ks_heapmake_uint64_t(size_t n, uint64_t *a);
void ks_heapup_128y(size_t n, ku128_t *a);
void ks_heapdown_128y(size_t i, size_t n, ku128_t *a);
void ks_heapmake_128y(size_t n, ku128_t *a);

void seq_char2nt6(int l, unsigned char *s);
void seq_reverse(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

uint64_t *fm6_seqsort(const rld_t *e, int n_threads);
int fm6_unitig(const struct __rld_t *e, int min_match, int n_threads, const uint64_t *sorted);
int fm6_ec_correct(const struct __rld_t *e, fmecopt_t *opt, const char *fn, int n_threads);
int fm6_remap(const char *fn, const rld_t *e, uint64_t *sorted, int skip, int min_pcv, int max_dist, int n_threads);
void mag_scaf_core(const rld_t *e, const char *fn, const fmscafopt_t *opt, int n_threads);

void fm_reverse_fmivec(fmintv_v *p);

uint64_t fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s, fmintv_t *k2, int *contained);

#endif
