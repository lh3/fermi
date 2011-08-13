#include "rld.h"
#include "fermi.h"

int sais(const unsigned char *T, int *SA, int n, int k);
int sais64(const unsigned char *T, int64_t *SA, int64_t n, int64_t k);

int fm_bwtgen(int asize, int64_t l, uint8_t *s)
{
	uint64_t i;
	if (l <= INT32_MAX) { // construct BWT for <2GB data
		int32_t *SA = malloc(l * 4);
		if (SA == 0) return -1;
		sais(s, SA, l, asize);
		for (i = 0; i < l; ++i) {
			if (SA[i] == 0) SA[i] = 0;
			else SA[i] = s[SA[i] - 1];
		}
		for (i = 0; i < l; ++i) s[i] = SA[i];
		free(SA);
	} else { // construct BWT for >=2GB data
		int64_t *SA = malloc(l * 8);
		if (SA == 0) return -1;
		sais64(s, SA, l, asize);
		for (i = 0; i < l; ++i) {
			if (SA[i] == 0) SA[i] = 0;
			else SA[i] = s[SA[i] - 1];
		}
		for (i = 0; i < l; ++i) s[i] = SA[i];
		free(SA);
	}
	return 0;
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
