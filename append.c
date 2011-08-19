#include <assert.h>
#include <stdio.h>
#include "utils.h"
#include "rld.h"

double cputime();
int sais(const unsigned char *T, int *SA, int n, int k);
void ks_introsort_uint64_t(size_t n, uint64_t a[]);

rld_t *fm_append(rld_t *e0, int len, const uint8_t *T)
{
	extern rld_t *fm_merge_from_SA(rld_t *e0, int len, const uint8_t *T, const int *SA, const uint64_t *rank_l);
	int c, k, *C, *p, *SA;
	uint64_t i, *oi, *rank_l;
	uint32_t *ws;
	rld_t *e;

	assert(T[len-1] == 0); // must be ended with a sentinel
	C = alloca(sizeof(int) * (e0->asize + 1));
	for (c = 0; c <= e0->asize; ++c) C[c] = 0;
	for (k = 0; k < len; ++k) ++C[T[k] + 1]; // marginal count
	for (c = 1; c <= e0->asize; ++c) C[c] += C[c-1]; // accumulative count
	ws = xmalloc((size_t)(len + 1) * 12);
	// set pointers
	rank_l = (uint64_t*)ws;
	SA = (int*)ws + 2 * (len + 1);
	// construct the suffix array
	sais(T, (int*)SA, len, e0->asize);
	// grab some memory from the stack
	p = alloca(sizeof(int) * e0->asize);
	oi = alloca(8 * e0->asize);
	// initialize the position array
	for (c = 0; c < e0->asize; ++c) p[c] = C[c+1] - 1; // point to the last element of each bucket
	// put the last sentinel 
	i = e0->mcnt[1] - 1;
	rank_l[p[0]--] = i;
	// compute the rank of long suffixes
	for (k = len - 2; k >= 0; --k) {
		if ((c = T[k]) != 0) {
			rld_rank1a(e0, i, oi);
			i = e0->cnt[c] + oi[c] - 1;
		} else i = e0->mcnt[1] - 1;
		rank_l[p[c]--] = i;
	}
	// sort the rank of long suffixes
	for (c = 1; c < e0->asize; ++c) // do not sort the sentinel bucket
		ks_introsort_uint64_t(C[c+1] - C[c], rank_l + C[c]);
	// merge to e0; e0 will be deallocated
	e = fm_merge_from_SA(e0, len, T, SA, rank_l);
	free(ws);
	return e;
}
