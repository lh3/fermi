#include <assert.h>
#include <stdlib.h>
#include "fermi.h"
#include "rld.h"

typedef struct {
	rlditr_t itr;
	int c;
	int64_t l;
} rlditr2_t;

uint32_t *fm_compute_gap0(const rld_t *e0, const rld_t *e1)
{
	uint32_t *gap;
	uint64_t k, l, *ok, *ol, i, j, x;
	int c = 0;
	ok = alloca(8 * e0->asize);
	ol = alloca(8 * e0->asize);
	gap = calloc(e0->mcnt[0], 4);
	x = e1->mcnt[1];
	k = l = --x; // get the last sentinel of e1
	j = i = e0->mcnt[1] - 1; // to modify gap[j]
	++gap[j];
	for (;;) {
		//printf("[%lld,%lld]@%d\n", k, l, c);
		rld_rank2a(e1, k - 1, l, ok, ol);
		for (c = 0; c < e1->asize; ++c)
			if (ok[c] < ol[c]) break;
		if (c == 0) {
			j = e0->mcnt[1] - 1;
			k = l = --x;
			if (x == (uint64_t)-1) break;
		} else {
			j = e0->cnt[c] + rld_rank11(e0, i, c) - 1;
			k = l = e1->cnt[c] + ok[c];
		}
		++gap[j];
		i = j;
	}
	return gap;
}

inline void rld_enc2(rld_t *e, rlditr2_t *itr2, int l, int c)
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
			l = rld_dec(e0, itr0, &c);
			assert(l); // the e0 stream should not be finished
			rld_enc2(e, itr, k < l? k : l, c);
		}
		*l0 = -k; *c0 = c;
	}
}

rld_t *fm_merge0(const rld_t *e0, const rld_t *e1)
{
	uint32_t *gap;
	uint64_t i, k;
	rlditr_t itr0, itr1;
	rlditr2_t itr;
	rld_t *e;
	int c0, c1;
	int64_t l0, l1;

	gap = fm_compute_gap0(e0, e1);
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr.itr, 0);
	itr.l = l0 = l1 = 0; itr.c = c0 = c1 = -1;
	rld_itr_init(e0, &itr0, 0);
	rld_itr_init(e1, &itr1, 0);
	for (i = k = 0; i < e0->mcnt[0]; ++i) {
		uint64_t g = gap[i];
		//printf("gap[%d]=%d\n", i, gap[i]);
		if (g) {
			//printf("%lld, %lld\n", g, k + 1);
			dec_enc(e, &itr, e0, &itr0, &l0, &c0, k + 1);
			dec_enc(e, &itr, e1, &itr1, &l1, &c1, g);
			k = 0;
		} else ++k;
	}
	if (k) dec_enc(e, &itr, e0, &itr0, &l0, &c0, k);
	assert(l0 == 0 && l1 == 0); // both e0 and e1 stream should be finished
	rld_enc(e, &itr.itr, itr.l, itr.c); // write the remaining symbols
	rld_enc_finish(e, &itr.itr);
	free(gap);
	return e;
}
