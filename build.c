#include <stdio.h>
#include <string.h>
#include "priv.h"

int fm_bwtgen(int asize, int64_t l, uint8_t *s)
{
	if (l <= INT32_MAX) return ksa_bwt(s, l, asize);
	else return ksa_bwt64(s, l, asize);
}

rld_t *fm_bwtenc(int asize, int sbits, int64_t l, const uint8_t *s)
{
	int c;
	int64_t i, k;
	rlditr_t itr;
	rld_t *e;

	e = rld_init(asize, sbits);
	rld_itr_init(e, &itr, 0);
	k = 1; c = s[0];
	for (i = 1; i < l; ++i) {
		if (s[i] != c) {
			rld_enc(e, &itr, k, c);
			c = s[i];
			k = 1;
		} else ++k;
	}
	rld_enc(e, &itr, k, c);
	rld_enc_finish(e, &itr);
	return e;
}

rld_t *fm_build(rld_t *e0, int asize, int sbits, int64_t l, uint8_t *s)
{
	rld_t *e;
	int64_t ori_l = e0? e0->mcnt[0] : 0;
	if (!e0) {
		fm_bwtgen(asize, l, s);
		e = fm_bwtenc(asize, sbits, l, s);
	} else e = fm_append(e0, l, s);
	if (fm_verbose >= 3) {
		int i;
		fprintf(stderr, "[M::%s] marginal counts: (", __func__);
		for (i = 0; i < e->asize1; ++i)
			fprintf(stderr, "%lld%s", (long long)e->mcnt[i], i == e->asize? ")" : ", ");
		fputc('\n', stderr);
	}
	assert(e->mcnt[0] == ori_l + l);
	return e;
}

rld_t *fm6_build2(int64_t l, const char *seq)
{
	int64_t i, j, beg;
	uint8_t *s;
	rld_t *e;
	s = calloc(l * 2, 1);
	for (i = j = beg = 0; i < l; ++i) { // generate the reverse complement
		s[j] = seq[i] < 6? seq[i] : seq_nt6_table[(int)seq[i]];
		if (s[j] == 0) {
			memcpy(s + j + 1, s + beg, j - beg);
			seq_revcomp6(j - beg, s + j + 1);
			j = beg = j - beg + 2 + j;
		} else ++j;
	}
	assert(j == l * 2);
	e = fm_build(0, 6, 3, l * 2, s);
	free(s);
	return e;
}
