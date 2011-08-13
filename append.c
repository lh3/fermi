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

static void gen_sufrank(const rld_t *e0, int len, const uint8_t *T, const int *C, sufrank_t *sr)
{
	uint64_t i, *oi;
	int c, k, *p;
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
		} else i = e0->mcnt[0] - 1;
		sr_rank(sr[p[c]]) = i;
		sr[p[c]--].suf = k;
	}
	for (c = 1; c < e0->asize; ++c) // do not sort the sentinel bucket
		ks_introsort(96, C[c+1] - C[c], sr + C[c]);
}

rld_t *fm_append(rld_t *e0, int len, uint8_t *T)
{
	int c, k, *C;
	rld_t *e;
	sufrank_t *sr;
	assert(T[len-1] == 0); // must be ended with a sentinel
	C = alloca(sizeof(int) * (e0->asize + 1));
	for (k = 0; k < len; ++k) ++C[T[k] + 1]; // marginal count
	for (c = 1; c <= e0->asize; ++c) C[c] += C[c-1]; // accumulative count
	sr = malloc(len * sizeof(sufrank_t));
	gen_sufrank(e0, len, T, C, sr);
	e = rld_init(e0->asize, e0->sbits);
	free(sr);
	return e;
}
