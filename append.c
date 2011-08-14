#include <assert.h>
#include "rld.h"
#include "ksort.h"

typedef struct {
	uint32_t rank[2]; // we do not use uint64_t because this will waste 25% of memory
	uint32_t suf;
} sufrank_t;

#define sr_rank(a) (*(uint64_t*)&(a))
#define sufrank_lt(a, b) (*(uint64_t*)&(a) < *(uint64_t*)&(b))
KSORT_INIT(96, sufrank_t, sufrank_lt)

int sais(const unsigned char *T, int *SA, int n, int k);
int sais_core(const unsigned char *T, int *SA, int fs, int n, int k, int cs);
void suffixsort(int *x, int *p, int n, int k, int l);
double cputime();

#include <stdio.h>

static void gen_sufrank(const rld_t *e0, int len, const uint8_t *T, const int *C, uint32_t *ws)
{
	sufrank_t *sr;
	uint64_t i, *oi, *rank_l, last;
	int c, k, l, *p;
	uint32_t *SA, *str;
	// set pointers
	sr = (sufrank_t*)ws;
	rank_l = (uint64_t*)ws; // rank_l and aux will replace sr at a later step
	SA = ws + 2 * (len + 1);
	str = SA + (len + 1);
	// allocate small memory
	p = alloca(sizeof(int) * e0->asize);
	oi = alloca(8 * e0->asize);
	// initialize the position array
	for (c = 0; c < e0->asize; ++c) p[c] = C[c+1] - 1; // point to the last element of each bucket
	// put the last sentinel 
	assert(T[len-1] == 0);
	i = e0->mcnt[1] - 1;
	sr_rank(sr[p[0]]) = i;
	sr[p[0]--].suf = len - 1;
	for (k = len - 2; k >= 0; --k) {
		if ((c = T[k]) != 0) {
			rld_rank1a(e0, i, oi);
			i = e0->cnt[c] + oi[c] - 1;
		} else i = e0->mcnt[1] - 1;
		sr_rank(sr[p[c]]) = i;
		sr[p[c]--].suf = k;
	}
	for (c = 1; c < e0->asize; ++c) // do not sort the sentinel bucket
		ks_introsort(96, C[c+1] - C[c], sr + C[c]);
	for (k = 0; k < C[1]; ++k) str[sr[k].suf] = k + 1;
	l = C[1]; last = sr_rank(sr[C[1]-1]);
	for (k = C[1]; k < len; ++k) {
		if (sr_rank(sr[k]) != last) last = sr_rank(sr[k]), ++l;
		str[sr[k].suf] = l;
	}
	for (k = 0; k < len; ++k) rank_l[k] = sr_rank(sr[k]);
	// The following lines gives another way to compute SA, but this takes more memory and we need to inverse SA to set str.
	double t = cputime();
	str[k] = 0; sais_core((unsigned char*)str, SA, 0, len + 1, l + 1, 4);
	fprintf(stderr, "%f\n", cputime() - t);
	t = cputime();
	suffixsort((int*)str, (int*)SA, len, l + 1, 0);
	fprintf(stderr, "%f\n", cputime() - t);
#if 1 // in principle, SA[] here should equal to the SA2[] computed directly from the input string
	int *SA2 = malloc(len * 4);
	t = cputime();
	sais(T, SA2, len, e0->asize);
	fprintf(stderr, "%f\n", cputime() - t);
	for (k = 0; k < len; ++k) assert(SA2[k] == SA[k+1]);
	free(SA2);
#endif
}

rld_t *fm_append(rld_t *e0, int len, uint8_t *T)
{
	int c, k, *C;
	rld_t *e;
	uint32_t *ws;
	assert(T[len-1] == 0); // must be ended with a sentinel
	C = alloca(sizeof(int) * (e0->asize + 1));
	for (c = 0; c <= e0->asize; ++c) C[c] = 0;
	for (k = 0; k < len; ++k) ++C[T[k] + 1]; // marginal count
	for (c = 1; c <= e0->asize; ++c) C[c] += C[c-1]; // accumulative count
	ws = malloc((len + 1) * 16);
	gen_sufrank(e0, len, T, C, ws);
	e = rld_init(e0->asize, e0->sbits);
	free(ws);
	return e;
}
