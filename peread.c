#include <stdio.h>
#include "fermi.h"
#include "utils.h"
#include "kvec.h"

#include "khash.h"
KHASH_DECLARE(64, uint64_t, uint64_t)

typedef khash_t(64) hash64_t;

void ks_introsort_128x(size_t n, fm128_t *a);

static hash64_t *build_hash(const fmnode_v *n)
{
	hash64_t *h;
	int64_t i, n_dropped = 0, n_dups = 0;
	int j, ret;
	khint_t k, l;

	h = kh_init(64);
	for (i = 0; i < n->n; ++i) {
		fmnode_t *p = &n->a[i];
		for (j = 0; j < p->mapping.n; ++j) {
			k = kh_put(64, h, p->mapping.a[j].x, &ret);
			kh_val(h, k) = ret == 0? (uint64_t)-1 : (i<<1|(p->mapping.a[j].y&1))^1;
		}
	}
	for (k = kh_begin(h); k != kh_end(h); ++k) {
		int to_drop = 0;
		if (!kh_exist(h, k)) continue;
		if (kh_val(h, k) != (uint64_t)-1) {
			l = kh_get(64, h, kh_key(h, k)^1);
			if (l == k) to_drop = 1;
			else if (l != kh_end(h)) {
				if (kh_val(h, l) == (uint64_t)-1) to_drop = 1;
			} else to_drop = 1;
		} else to_drop = 1, ++n_dups;
		if (to_drop) {
			kh_del(64, h, k);
			++n_dropped;
		}
	}
	if (fm_verbose >= 3)
		fprintf(stderr, "[%s] dropped %ld reads, including %ld duplicates; %d reads remain\n",
				__func__, (long)n_dropped, (long)n_dups, kh_size(h));
	return h;
}

static void collect_pairs(const fmnode_v *n, const hash64_t *h, fm128_v *pairs)
{
	size_t i, m;
	int j;
	khint_t k;
	fm128_t z, *r;
	pairs->n = 0;
	for (i = 0; i < n->n; ++i) {
		fmnode_t *p = &n->a[i];
		for (j = 0; j < p->mapping.n; ++j) {
			k = kh_get(64, h, p->mapping.a[j].x);
			if (k == kh_end(h)) continue;
			z.x = kh_val(h, k)<<8;
			k = kh_get(64, h, p->mapping.a[j].x^1);
			z.y = kh_val(h, k);
			if (z.x>>8 < z.y) kv_push(fm128_t, *pairs, z);
		}
	}
	ks_introsort_128x(pairs->n, pairs->a);
	r = &pairs->a[0];
	for (i = 1; i < pairs->n; ++i) {
		fm128_t *q = &pairs->a[i];
		if (q->x>>8 == r->x>>8 && q->y == r->y) {
			if ((r->x&0xff) != 0xff) ++r->x;
			q->x = q->y = (uint64_t)-1;
		} else r = q, ++r->x;
	}
	for (i = m = 0; i < pairs->n; ++i)
		if (pairs->a[i].x != (uint64_t)-1 && (pairs->a[i].x&0xff) > 1)
			pairs->a[m++] = pairs->a[i];
	pairs->n = m;
	for (i = 0; i < pairs->n; ++i) {
		printf("%lld, %lld, %lld\n", pairs->a[i].x>>8, pairs->a[i].y, pairs->a[i].x&0xff);
	}
}

int msg_peread(const msg_t *g, int max_dist)
{
	hash64_t *h;
	fm128_v pairs;
	kv_init(pairs);
	h = build_hash(&g->nodes);
	collect_pairs(&g->nodes, h, &pairs);
	kh_destroy(64, h);
	free(pairs.a);
	return 0;
}
