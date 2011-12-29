#ifndef FM_PRIV_H
#define FM_PRIV_H

#include "rld.h"
#include "fermi.h"
#include "utils.h"

int ksa_sa(const unsigned char *T, int *SA, int n, int k);
int ksa_bwt(unsigned char *T, int n, int k);
int ksa_bwt64(unsigned char *T, int64_t n, int k);

void seq_char2nt6(int l, unsigned char *s);
void seq_reverse(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

uint64_t *fm6_seqsort(const rld_t *e, int n_threads);
int fm6_unitig(const struct __rld_t *e, int min_match, int n_threads, const uint64_t *sorted);
int fm6_ec_correct(const struct __rld_t *e, const fmecopt_t *opt, const char *fn, int n_threads);

void fm_reverse_fmivec(fmintv_v *p);

uint64_t fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s, fmintv_t *k2, int *contained);

void msg_write_node(const fmnode_t *p, long id, kstring_t *out);
void msg_nodecpy(fmnode_t *dst, const fmnode_t *src);
void msg_join_unambi(msg_t *g);

#endif
