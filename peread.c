#include <stdio.h>
#include <limits.h>
#include "fermi.h"
#include "utils.h"
#include "kvec.h"

#include "khash.h"
KHASH_DECLARE(64, uint64_t, uint64_t)

typedef khash_t(64) hash64_t;

void ks_introsort_128x(size_t n, fm128_t *a);
void ks_heapup_128y(size_t n, fm128_t *a);
void ks_heapdown_128y(size_t i, size_t n, fm128_t *a);

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
			k = kh_put(64, h, p->mapping.a[j].x, &ret); // the key is the read index
			if (ret) { // absent from the hash table
				uint64_t tmp = p->mapping.a[j].y;
				uint32_t pos = (tmp&1)? tmp<<32>>33 : tmp>>32; // pos of the end of a fragment
				kh_val(h, k) = i<<32 | ((tmp&1)^1)<<31 | pos;
			} else kh_val(h, k) = (uint64_t)-1; // a read with multiple occurrences; drop it
		}
	}
	for (k = kh_begin(h); k != kh_end(h); ++k) { // exclude unpaired reads
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

static void collect_pairs(fmnode_v *n, fm128_v *pairs) // n->a[].aux to be modified
{
	size_t i;
	fm128_t z;
	hash64_t *h;
	h = build_hash(n);
	pairs->n = 0;
	for (i = 0; i < n->n; ++i) {
		int j;
		fmnode_t *p = &n->a[i];
		p->aux[0] = p->aux[1] = INT_MAX;
		for (j = 0; j < p->mapping.n; ++j) {
			khint_t k;
			k = kh_get(64, h, p->mapping.a[j].x);
			if (k == kh_end(h)) continue;
			z.x = kh_val(h, k)>>31<<32;
			z.y = kh_val(h, k)<<33>>33<<32;
			k = kh_get(64, h, p->mapping.a[j].x^1);
			z.x |= kh_val(h, k)>>31;
			z.y |= kh_val(h, k)<<33>>33;
			if (z.x>>32 < z.x<<32>>32) kv_push(fm128_t, *pairs, z);
		}
		free(p->mapping.a);
		p->mapping.a = 0;
		p->mapping.n = p->mapping.m = 0;
	}
	kh_destroy(64, h);
	ks_introsort_128x(pairs->n, pairs->a);
}

static void index_pairs(const fm128_v *pairs, fm64_v *pidx)
{
	uint32_t i, beg;
	for (i = 1, beg = 0; i < pairs->n; ++i)
		if (pairs->a[i].x != pairs->a[beg].x) {
			kv_push(uint64_t, *pidx, (uint64_t)beg<<32 | i);
			beg = i;
		}
	kv_push(uint64_t, *pidx, (uint64_t)beg<<32 | i);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] collected %ld read pairs and %ld unitig pairs\n", __func__, (long)pairs->n, (long)pidx->n);
}

static inline uint64_t get_idd(hash64_t *h, uint64_t k)
{
	khint_t iter;
	iter = kh_get(64, h, k);
	return iter == kh_end(h)? (uint64_t)-1 : kh_val(h, iter);
}

typedef struct {
	fm128_v heap, stack, rst;
	fm64_v walk;
	int is_multi;
} aux_t;

//static int walk(msg_t *g, const hash64_t *h, size_t idd[2], int max_dist, aux_t *a)
static int walk(msg_t *g, size_t idd[2], int max_dist, aux_t *a)
{ // FIXME: the algorithm can be improved but will be more complicated.
	fm128_t *q;
	fm128_v *r;
	fmnode_t *p, *w;
	int i, n_nei[2];
	uint64_t end, start;

	// stack -- .x: id+direction (idd); .y: parent
	// heap  -- .x: position in stack;  .y: distance from the end of start (can be negative)

	// initialize
	a->heap.n = a->stack.n = a->rst.n = a->walk.n = a->is_multi = 0;
	for (i = 0; i < 2; ++i) {
		p = &g->nodes.a[i];
		n_nei[i] = p->nei[idd[i]&1].n;
	}
	if (n_nei[0] < n_nei[1]) start = idd[0], end = idd[1];
	else start = idd[1], end = idd[0];
	kv_pushp(fm128_t, a->stack, &q);
	q->x = start, q->y = (uint64_t)-1;
	kv_pushp(fm128_t, a->heap, &q);
	q->x = 0, q->y = -g->nodes.a[start>>1].l; // note that "128y" compares int64_t instead of uint64_t
	g->nodes.a[start>>1].aux[0] = 1;
	// shortest path
	while (a->heap.n) {
		// pop up the best node
		fm128_t z = a->heap.a[0];
		a->heap.a[0] = a->heap.a[--a->heap.n];
		ks_heapdown_128y(0, a->heap.n, a->heap.a);
		// push to the heap
		p = &g->nodes.a[a->stack.a[z.x].x>>1];
		r = &p->nei[a->stack.a[z.x].x&1];
		for (i = 0; i < r->n; ++i) {
			uint64_t u = get_idd(g->h, r->a[i].x);
			int64_t dist = (int64_t)z.y + p->l - (int64_t)r->a[i].y;
			w = &g->nodes.a[u>>1];
			if (dist < max_dist) {
				if (w->aux[0] != INT_MAX) { // visited before
					++w->aux[0];
					if (fm_verbose >= 100) printf("multi: [%lld,%lld]->[%lld,%lld]\n", p->k[0], p->k[1], w->k[0], w->k[1]);
					if (a->rst.n) break;
				} else {
					kv_pushp(fm128_t, a->heap, &q);
					q->x = a->stack.n, q->y = dist;
					w->aux[0] = 1;
					ks_heapup_128y(a->heap.n, a->heap.a);
					kv_pushp(fm128_t, a->stack, &q);
					q->x = u^1, q->y = z.x;
					if (fm_verbose >= 100) printf("[%lld,%lld]->[%lld,%lld]\t%lld\n", p->k[0], p->k[1], w->k[0], w->k[1], dist);
					if (u == end) { // reach the end
						kv_pushp(fm128_t, a->rst, &q);
						q->x = a->stack.n - 1, q->y = dist;
					}
				}
			}
		}
		if (i != r->n) break; // found 3 paths
	}
	if (a->rst.n == 0) { // no path
		for (i = 0; i < a->stack.n; ++i)
			p = &g->nodes.a[a->stack.a[i].x>>1], p->aux[0] = p->aux[1] = INT_MAX;
		return INT_MIN;
	}
	// backtrace
	end = a->rst.a[0].x;
	do {
		if (g->nodes.a[a->stack.a[end].x>>1].aux[0] > 1) a->is_multi = 1;
		kv_push(uint64_t, a->walk, a->stack.a[end].x);
		end = a->stack.a[end].y;
	} while (end != (uint64_t)-1);
	for (i = 0; i < a->walk.n>>1; ++i) // reverse
		end = a->walk.a[i], a->walk.a[i] = a->walk.a[a->walk.n - 1 - i], a->walk.a[a->walk.n - 1 - i] = end;
	for (i = 0; i < a->stack.n; ++i)
		p = &g->nodes.a[a->stack.a[i].x>>1], p->aux[0] = p->aux[1] = INT_MAX;
	return (int)((int64_t)a->rst.a[0].y);
}

int msg_peread(msg_t *g, double avg, double std)
{
	int64_t i;
	int j, max_dist = (int)(avg + std * 2 + .499);
	fm128_v pairs;
	fm64_v pidx;
	aux_t a;

	kv_init(pairs); kv_init(pidx);
	memset(&a, 0, sizeof(aux_t));

	collect_pairs(&g->nodes, &pairs);
	index_pairs(&pairs, &pidx);
	for (i = 0; i < pidx.n; ++i) {
		int dist;
		size_t idd[2];
		fm128_t *q;
		q = &pairs.a[pidx.a[i]>>32];
		idd[0] = q->x>>32; idd[1] = q->x<<32>>32;
		fm_verbose = (idd[0] == 144 && idd[1] == 268)? 1000 : 3;
		dist = walk(g, idd, max_dist, &a);
		printf("***\t%d\t%lld[%lld]\t%lld[%lld]\t", (int)((pidx.a[i]<<32>>32) - (pidx.a[i]>>32)), idd[0], g->nodes.a[idd[0]>>1].k[idd[0]&1], idd[1], g->nodes.a[idd[1]>>1].k[idd[1]&1]);
		if (dist == INT_MIN) {
			printf("none");
		} else {
			printf("%d\t%d\t", a.is_multi, dist);
			for (j = 0; j < a.walk.n; ++j) {
				if (j) putchar(',');
				printf("%lld", a.walk.a[j]);
			}
		}
		putchar('\n');
		fflush(stdout);
	}

	free(a.walk.a); free(a.rst.a); free(a.stack.a); free(a.heap.a);
	free(pairs.a); free(pidx.a);
	return 0;
}
