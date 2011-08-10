#include <assert.h>
#include <stdlib.h>
#include "fermi.h"
#include "rld.h"

typedef struct {
	rlditr_t itr;
	int l, c;
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
	while (x != (uint64_t)-1 || c) {
		//printf("[%lld,%lld]@%d, %lld\n", k, l, c, x);
		rld_rank2a(e1, k - 1, l, ok, ol);
		for (c = 0; c < e1->asize; ++c)
			if (ok[c] < ol[c]) break;
		if (c == 0) {
			j = e0->mcnt[1] - 1;
			k = l = --x;
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
	if (itr2->c != c) {
		if (itr2->l) rld_enc(e, &itr2->itr, itr2->l, itr2->c);
		itr2->l = l; itr2->c = c;
	} else itr2->l += l;
}

rld_t *fm_merge0(const rld_t *e0, const rld_t *e1)
{
	uint32_t *gap;
	uint64_t i, k;
	rlditr_t itr0, itr1;
	rlditr2_t itr;
	rld_t *e;
	int l0, c0, l1, c1, l, c;

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
			if (l0 >= (int)k) {
				rld_enc2(e, &itr, k, c0);
				l0 -= k;
			} else {
				rld_enc2(e, &itr, l0, c0);
				k -= l0;
				for (;;) {
					l = rld_dec(e0, &itr0, &c);
					if (l >= (int)k) {
						rld_enc2(e, &itr, k, c);
						l0 = l - k; c0 = c;
						break;
					} else rld_enc2(e, &itr, l, c);
					k -= l;
				}
			}
			k = g;
			if (l1 >= (int)k) {
				rld_enc2(e, &itr, k, c1);
				l1 -= k;
			} else {
				rld_enc2(e, &itr, l1, c1);
				k -= l1;
				for (;;) {
					l = rld_dec(e1, &itr1, &c);
					if (l >= (int)k) {
						rld_enc2(e, &itr, k, c);
						l1 = l - k; c1 = c;
						break;
					} else rld_enc2(e, &itr, l, c);
					k -= l;
				}
			}
		} else ++k;
	}
	// FIXME: deal with the remaining l0 and l1
	rld_enc(e, &itr.itr, itr.l, itr.c);
	rld_enc_finish(e, &itr.itr);
	free(gap);
	return e;
}
