#include "priv.h"
#include "mog.h"
#include "kvec.h"

#define edge_mark_del(_x) ((_x).x = (uint64_t)-2, (_x).y = 0)
#define edge_is_del(_x)   ((_x).x == (uint64_t)-2 || (_x).y == 0)

typedef struct {
	uint32_t id;
	int8_t finished[2];
	int cnt[2];
	int32_t d[2], d2[2];
	uint32_t p[2], p2[2];
} trinfo_t;

static trinfo_t g_ti_null = { UINT32_MAX, {0, 0}, {0, 0}, {INT32_MAX, INT32_MAX}, {INT32_MAX, INT32_MAX}, {UINT32_MAX, UINT32_MAX}, {UINT32_MAX, UINT32_MAX}};

typedef struct {
	int n, m;
	trinfo_t **buf;
} tipool_t;

typedef struct {
	tipool_t pool;
	ku64_v heap, stack;
} bblaux_t;

#define tiptr(p) ((trinfo_t*)(p)->ptr)

static inline trinfo_t *tip_alloc(tipool_t *pool, uint32_t id)
{
	trinfo_t *p;
	if (pool->n == pool->m) {
		int i, new_m = pool->m? pool->m<<1 : 256;
		pool->buf = realloc(pool->buf, new_m);
		for (i = pool->m; i < new_m; ++i)
			pool->buf[i] = malloc(sizeof(trinfo_t));
	}
	p = pool->buf[pool->n++];
	*p = g_ti_null;
	p->id = id;
	return p;
}

void mog_vh_bbl_detect(mog_t *g, uint64_t idd, int max_vtx, int max_dist, bblaux_t *a)
{
	int i, shift;
	int64_t term = -1;
	mogv_t *p, *q;

	a->heap.n = a->stack.n = a->pool.n = 0;
	p = &g->v.a[idd>>1];
	if (p->len < 0 || p->nei[idd&1].n < 2) return; // stop if p is deleted or it has 0 or 1 neighbor
	shift = p->len; // shift to make sure the distance in the higher 32-bit is always non-negative
	p->ptr = tip_alloc(&a->pool, idd>>1);
	tiptr(p)->d[(idd&1)^1] = -shift;
	kv_push(uint64_t, a->heap, idd^1);
	// Modified Dijkstra's algorithm for top 2 shortest paths
	while (a->heap.n) {
		uint64_t x, w;
		ku128_v *r;
		// get the closest vertex
		x = a->heap.a[0]; a->heap.a[0] = a->heap.a[--a->heap.n];
		ks_heapdown_uint64_t(0, a->heap.n, a->heap.a);
		p = &g->v.a[((uint32_t)x)>>1];
		r = &p->nei[(x&1)^1]; // we will look the the nighbors from the other end of the unitig
		if (tiptr(p)->finished[x&1]) continue;
		// test if stop
		if ((int)(x>>32) - shift > max_dist || r->n == 0) { // we come to a dead end
			if (term < 0) term = (uint32_t)x;
			else { term = -2; break; } // multiple dead ends
		}
		// set the distance to p's neighbors
		for (i = 0; i < r->n; ++i) {
			int d, d2;
			if (edge_is_del(r->a[i])) continue;
			w = mog_tid2idd(g->h, r->a[i].x);
			q = &g->v.a[w>>1];
			if (q->ptr == 0) q->ptr = tip_alloc(&a->pool, w>>1);
			++tiptr(q)->cnt[w&1];
			d = tiptr(p)->d[x&1] + p->len - r->a[i].y; // distance to the end of q
			// test and possibly update the tentative distance
			if (d < tiptr(q)->d[w&1]) {
				// move the best to the 2nd best
				tiptr(q)->d2[w&1] = tiptr(q)->d[w&1];
				tiptr(q)->p2[w&1] = tiptr(q)->p[w&1];
				// update the bes
				tiptr(q)->d[w&1] = d;
				tiptr(q)->p[w&1] = (uint32_t)x;
				// push to the heap
				kv_push(uint64_t, a->heap, (uint64_t)(d + shift) << 32 | w);
				ks_heapup_uint64_t(a->heap.n, a->heap.a);
				// test and update d2
				d2 = tiptr(p)->d2[x&1] + p->len - r->a[i].y;
				if (d2 < tiptr(q)->d2[w&1]) {
					tiptr(q)->d2[w&1] = d2;
					tiptr(q)->p2[w&1] = (uint32_t)x;
				}
			} else if (d < tiptr(q)->d2[w&1]) {
				tiptr(q)->d2[w&1] = d;
				tiptr(q)->p2[w&1] = (uint32_t)x;
			}
		}
		tiptr(p)->finished[x&1] = 1;
	}
	// test if the bubble is complete
}
