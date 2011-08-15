#include <assert.h>
#include <stdio.h>
#include "rld.h"
#include "ksort.h"

typedef struct {
	uint32_t rank[2]; // we do not use uint64_t because this will waste 25% of memory
	uint32_t suf;
} sufrank_t;

double cputime();
void ks_introsort_uint64_t(size_t, uint64_t*); // defined in merge.c
int sais(const unsigned char *T, int *SA, int n, int k);

#define sr_rank(a) (*(uint64_t*)&(a))
#define sufrank_lt(a, b) (*(uint64_t*)&(a) < *(uint64_t*)&(b))
KSORT_INIT(96, sufrank_t, sufrank_lt)

void qsufsort_mod(int *V, int *I, int numChar, int largestInputSymbol, int smallestInputSymbol);

rld_t *fm_append(rld_t *e0, int len, const uint8_t *T)
{
	int c, k, *C, *p, *rank_s, *aux;
	uint64_t i, *oi, *rank_l;
	uint32_t *ws;
	rld_t *e;
	sufrank_t *sr;

	assert(T[len-1] == 0); // must be ended with a sentinel
	C = alloca(sizeof(int) * (e0->asize + 1));
	for (c = 0; c <= e0->asize; ++c) C[c] = 0;
	for (k = 0; k < len; ++k) ++C[T[k] + 1]; // marginal count
	for (c = 1; c <= e0->asize; ++c) C[c] += C[c-1]; // accumulative count
	ws = malloc((len + 1) * 16);
	// set pointers
	sr = (sufrank_t*)ws;
	rank_l = (uint64_t*)ws;
	rank_s = (int*)ws + 2 * (len + 1);
	aux = rank_s + len + 1;
	// grab some memory from the stack
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
	// sort the rank of long suffixes
	for (c = 1; c < e0->asize; ++c) // do not sort the sentinel bucket
		ks_introsort(96, C[c+1] - C[c], sr + C[c]);
	for (k = 0; k < len; ++k) aux[k] = sr[k].suf, rank_l[k] = sr_rank(sr[k]);
	{
		int i, li, s, r, f;
		uint64_t lr;
		li = len - 1;
		lr = rank_l[li];
		s = aux[li];
		rank_s[s] = li + 1;
		f = C[c = e0->asize - 1];
		for (i = len - 1; i--; ) {
			r = rank_l[i], s = aux[i];
			if (i < f) {
				if (li >= f) ++lr;
				f = C[--c];
			}
			if (r != lr || c == 0) {
				li = i, lr = r;
				rank_s[s] = i + 1;
			} else rank_s[s] = li + 1;
		}
		rank_s[len] = 0;
		aux[len] = aux[len - 1]; aux[0] = -1;
		for (i = len - 2; i >= C[1]; --i)
			aux[i+1] = (rank_l[i] != rank_l[i+1] && rank_l[i] != rank_l[i-1])? -1 : aux[i];
		for (; i >= 0; --i) aux[i+1] = -1;
		qsufsort_mod(rank_s, aux, len, len, 1);
		//suffixsort(rank_s, aux, len, len+1, 1);
#if 0
		int *SA = malloc((len+1) * sizeof(int));
		sais(T, SA, len, e0->asize);
		for (i = 0; i < len; ++i) assert(i == SA[rank_s[i]-1]);
#endif
	}
	free(ws);
	exit(1);
	return e;
}

rld_t *fm_appenda(rld_t *e0, int len, const uint8_t *T)
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
	ws = malloc((len + 1) * 12);
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
		ks_introsort(uint64_t, C[c+1] - C[c], rank_l + C[c]);
	// merge to e0; e0 will be deallocated
	e = fm_merge_from_SA(e0, len, T, SA, rank_l);
	free(ws);
	return e;
}
