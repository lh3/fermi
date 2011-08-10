#include <assert.h>
#include <stdlib.h>
#include "fermi.h"
#include "rld.h"

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
		printf("[%lld,%lld]@%d, %lld\n", k, l, c, x);
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

rld_t *fm_merge0(const rld_t *e0, const rld_t *e1)
{
	uint32_t *gap;
	int i;
	gap = fm_compute_gap0(e0, e1);
	for (i = 0; i < e0->mcnt[0]; ++i) printf("gap[%d]=%d\n", i, gap[i]);
	free(gap);
	return 0;
}
