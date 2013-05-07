#include "fermi.h"
#include "rld.h"

#include "klist.h"
#define myfree(x)
KLIST_INIT(intv, ku128_t, myfree)

void fm6_bfs_traverse(const rld_t *e, uint64_t k0, uint64_t l0)
{
	uint64_t *bl, *br;
	klist_t(intv) *list;
	ku128_t *p;

	list = kl_init(intv);
	bl = calloc((e->mcnt[0] + 64) / 64, 8);
	br = calloc((e->mcnt[0] + 64) / 64, 8);
	p = kl_pushp(intv, list);
	p->x = k0; p->y = l0;
	while (list->size) {
		uint64_t yl, yr, *pl, *pr, zl, zr;
		uint64_t ok[6], ol[6];
		ku128_t x;
		int c;

		kl_shift(intv, list, &x);
		printf("[%ld,%ld)\n", (long)x.x, (long)x.y);
		rld_rank2a(e, x.x - 1, x.y - 1, ok, ol);
		for (c = 1; c < e->asize; ++c) {
			if (ok[c] == ol[c]) continue;
			x.x = e->cnt[c] + ok[c];
			x.y = e->cnt[c] + ol[c];
			pl = bl + (x.x>>6); zl = 1ULL << (x.x&0x3f);
			pr = br + (x.y>>6); zr = 1ULL << (x.y&0x3f);
			yl = __sync_fetch_and_or(pl, zl);
			yr = __sync_fetch_and_or(pr, zr);
			if ((yl&zl) && (yr&zr)) continue;
			p = kl_pushp(intv, list);
			p->x = e->cnt[c] + ok[c];
			p->y = e->cnt[c] + ol[c];
		}
	}
	free(bl); free(br);
	kl_destroy(intv, list);
}
