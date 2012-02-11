#include <limits.h>
#include "priv.h"
#include "mog.h"
#include "kvec.h"

#define edge_mark_del(_x) ((_x).x = (uint64_t)-2, (_x).y = 0)
#define edge_is_del(_x)   ((_x).x == (uint64_t)-2 || (_x).y == 0)

typedef struct {
	uint32_t id;
	int cnt[2];
	int n[2][2], d[2][2];
	uint32_t v[2][2];
} trinfo_t;

typedef struct {
	int n, m;
	trinfo_t **buf;
} tipool_t;

typedef struct {
	tipool_t pool;
	ku64_v stack;
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
	memset(p, 0, sizeof(trinfo_t));
	p->id = id;
	return p;
}

void mog_vh_bbl_detect(mog_t *g, uint64_t idd, int max_vtx, int max_dist, bblaux_t *a)
{
	int i, n_pending = 0;
	mogv_t *p, *q;

	a->stack.n = a->pool.n = 0;
	p = &g->v.a[idd>>1];
	if (p->len < 0 || p->nei[idd&1].n < 2) return; // stop if p is deleted or it has 0 or 1 neighbor
	p->ptr = tip_alloc(&a->pool, idd>>1);
	tiptr(p)->d[(idd&1)^1][0] = tiptr(p)->d[(idd&1)^1][1] = -p->len;
	kv_push(uint64_t, a->stack, idd^1);
	while (a->stack.n) {
		uint64_t x, y;
		ku128_v *r;
		if (a->stack.n == 1 && a->stack.a[0] != (idd^1) && n_pending == 0) break; // found the other end of the bubble
		x = kv_pop(a->stack);
		p = &g->v.a[x>>1];
		r = &p->nei[(x&1)^1]; // we will look the the neighbors from the other end of the unitig
		if (tiptr(p)->d[x&1][0] > max_dist || tiptr(p)->d[x&1][1] > max_dist || r->n == 0) break; // we failed
		// set the distance to p's neighbors
		for (i = 0; i < r->n; ++i) {
			int nsr;
			if (edge_is_del(r->a[i])) continue;
			y = mog_tid2idd(g->h, r->a[i].x);
			q = &g->v.a[y>>1];
			if (q->ptr == 0) { // has not been attempted
				q->ptr = tip_alloc(&a->pool, y>>1), ++n_pending;
				mog_v128_clean(&q->nei[y&1]); // make sure there are no deleted edges
			}
			nsr = tiptr(p)->n[x&1][0] + q->nsr;
			// test and possibly update the tentative distance
			if (nsr > tiptr(q)->n[y&1][0]) { // then move the best to the 2nd best and update the best
				tiptr(q)->n[y&1][1] = tiptr(q)->n[y&1][0]; tiptr(q)->n[y&1][0] = nsr;
				tiptr(q)->v[y&1][1] = tiptr(q)->v[y&1][0]; tiptr(q)->v[y&1][0] = x;
				tiptr(q)->d[y&1][1] = tiptr(q)->d[y&1][0]; tiptr(q)->d[y&1][0] = tiptr(p)->d[x&1][0] + p->len - r->a[i].y;
				nsr = tiptr(p)->n[x&1][1] + q->nsr; // now nsr is the 2nd best
			}
			if (nsr > tiptr(q)->n[y&1][1]) { // update the 2nd best
				tiptr(q)->n[y&1][1] = nsr, tiptr(q)->v[y&1][1] = x;
				tiptr(q)->d[y&1][1] = tiptr(p)->d[x&1][1] + p->len - r->a[i].y;
			}
			if (++tiptr(q)->cnt[y&1] == q->nei[y&1].n) { // all q's predecessors have been processed; then push
				kv_push(uint64_t, a->stack, y);
				--n_pending;
			}
		}
	}
	for (i = 0; i < a->pool.n; ++i) // reset p->ptr
		g->v.a[a->pool.buf[i]->id].ptr = 0;
}
