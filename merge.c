#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "fermi.h"
#include "rld.h"
#include "utils.h"

double cputime();

uint64_t *fm_compute_gap_bits(const rld_t *e0, const rld_t *e1, int n_threads);

typedef struct {
	rlditr_t itr;
	int c;
	int64_t l;
} rlditr2_t;

// This is a clever version of rld_enc(): adjacent runs are guaranteed to be different.
inline void rld_enc2(rld_t *e, rlditr2_t *itr2, int64_t l, int c)
{
	if (l == 0) return;
	if (itr2->c != c) {
		if (itr2->l) rld_enc(e, &itr2->itr, itr2->l, itr2->c);
		itr2->l = l; itr2->c = c;
	} else itr2->l += l;
}
// take k symbols from e0 and write it to e; l0 number of c0 symbols are pending before writing
static inline void dec_enc(rld_t *e, rlditr2_t *itr, const rld_t *e0, rlditr_t *itr0, int64_t *l0, int *c0, int64_t k)
{
	if (*l0 >= k) { // there are more pending symbols
		rld_enc2(e, itr, k, *c0);
		*l0 -= k; // l0-k symbols remains
	} else { // use up all pending symbols
		int c = -1; // to please gcc
		int64_t l;
		rld_enc2(e, itr, *l0, *c0); // write all pending symbols
		k -= *l0;
		for (; k > 0; k -= l) { // we always go into this loop because l0<k
			l = rld_dec(e0, itr0, &c, 1);
			assert(l); // the e0 stream should not be finished
			rld_enc2(e, itr, k < l? k : l, c);
		}
		*l0 = -k; *c0 = c;
	}
}

rld_t *fm_merge(rld_t *e0, rld_t *e1, int use_hash, int n_threads)
{
	uint64_t i, n = e0->mcnt[0] + e1->mcnt[0], *bits;
	int64_t l0, l1;
	int c0, c1;
	rlditr_t itr0, itr1;
	rlditr2_t itr;
	rld_t *e;

	// compute the gap array
	bits = fm_compute_gap_bits(e0, e1, n_threads);
	free(e0->frame); free(e1->frame); // deallocate the rank indexes of e0 and e1; they are not needed any more
	e0->frame = e1->frame = 0;
	// initialize the FM-index to be returned, and all the three iterators
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr.itr, 0);
	itr.l = l0 = l1 = 0; itr.c = c0 = c1 = -1;
	rld_itr_init(e0, &itr0, 0);
	rld_itr_init(e1, &itr1, 0);
	{
		uint64_t k = 1;
		int last = bits[0]&1;
		for (i = 1; i < n; ++i) {
			int c = bits[i>>6]>>(i&0x3f)&1;
			if (c != last) {
				if (last == 0) dec_enc(e, &itr, e0, &itr0, &l0, &c0, k);
				else dec_enc(e, &itr, e1, &itr1, &l1, &c1, k);
				last = c; k = 1;
			} else ++k;
		}
		if (k) {
			if (last == 0) dec_enc(e, &itr, e0, &itr0, &l0, &c0, k);
			else dec_enc(e, &itr, e1, &itr1, &l1, &c1, k);
		}
		free(bits);
	}
	// finalize the merge
	assert(l0 == 0 && l1 == 0); // both e0 and e1 stream should be finished
	rld_enc(e, &itr.itr, itr.l, itr.c); // write the remaining symbols in the iterator
	rld_destroy(e0); rld_destroy(e1);
	rld_enc_finish(e, &itr.itr);
	return e;
}

rld_t *fm_merge_from_SA(rld_t *e0, int len, const uint8_t *T, const int *SA, const int64_t *rank_l)
{
	int64_t l0, last = -1;
	int c0, i;
	rlditr_t itr0;
	rlditr2_t itr;
	rld_t *e;

	free(e0->frame); e0->frame = 0;
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr.itr, 0);
	itr.l = l0 = 0; itr.c = c0 = -1;
	rld_itr_init(e0, &itr0, 0);
	for (i = 0; i < len; ++i) {
		if (rank_l[i] != last) {
			dec_enc(e, &itr, e0, &itr0, &l0, &c0, rank_l[i] - last);
			last = rank_l[i];
		}
		rld_enc2(e, &itr, 1, SA[i]? T[SA[i]-1] : 0);
	}
	if (last != e0->mcnt[0] - 1)
		dec_enc(e, &itr, e0, &itr0, &l0, &c0, e0->mcnt[0] - 1 - last);
	rld_enc(e, &itr.itr, itr.l, itr.c); // write the remaining symbols in the iterator
	rld_destroy(e0);
	rld_enc_finish(e, &itr.itr);
	return e;
}
