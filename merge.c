#include <assert.h>
#include <stdlib.h>
#include "fermi.h"
#include "rld.h"
#include "kvec.h"

typedef struct {
	rlditr_t itr;
	int c;
	int64_t l;
} rlditr2_t;

typedef struct {
	int32_t *gap;
	kvec_t(int64_t) a;
} gaparr_t;

#define GAP_MAX INT32_MAX

static gaparr_t *compute_gap_array(const rld_t *e0, const rld_t *e1)
{
	gaparr_t *g;
	uint64_t k, l, *ok, *ol, i, j, x;
	int c = 0;
	g = calloc(1, sizeof(gaparr_t));
	ok = alloca(8 * e0->asize);
	ol = alloca(8 * e0->asize);
	g->gap = calloc(e0->mcnt[0], 4);
	x = e1->mcnt[1];
	k = l = --x; // get the last sentinel of e1
	j = i = e0->mcnt[1] - 1; // to modify gap[j]
	++g->gap[j];
	for (;;) {
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
		if (g->gap[j] < 0) {
			++g->a.a[-g->gap[j] - 1];
		} else if (g->gap[j] == GAP_MAX) {
			kv_push(int64_t, g->a, 1ll + GAP_MAX);
			g->gap[j] = -(int64_t)g->a.n;
		} else ++g->gap[j];
		i = j;
	}
	return g;
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

rld_t *fm_merge_array(rld_t *e0, rld_t *e1, const char *fn)
{
	gaparr_t *gap;
	uint64_t i, k;
	rlditr_t itr0, itr1;
	rlditr2_t itr;
	rld_t *e;
	int c0, c1;
	int64_t l0, l1;

	gap = compute_gap_array(e0, e1);
	e = rld_init(e0->asize, e0->sbits, fn);
	rld_itr_init(e, &itr.itr, 0);
	itr.l = l0 = l1 = 0; itr.c = c0 = c1 = -1;
	rld_itr_init(e0, &itr0, 0);
	rld_itr_init(e1, &itr1, 0);
	for (i = k = 0; i < e0->mcnt[0]; ++i) {
		int64_t g = gap->gap[i] < 0? gap->a.a[-gap->gap[i] - 1] : gap->gap[i];
		if (g) {
			//printf("gap[%lld]=%lld\n", i, g);
			dec_enc(e, &itr, e0, &itr0, &l0, &c0, k + 1);
			dec_enc(e, &itr, e1, &itr1, &l1, &c1, g);
			k = 0;
		} else ++k;
	}
	if (k) dec_enc(e, &itr, e0, &itr0, &l0, &c0, k);
	assert(l0 == 0 && l1 == 0); // both e0 and e1 stream should be finished
	rld_enc(e, &itr.itr, itr.l, itr.c); // write the remaining symbols
	rld_enc_finish(e, &itr.itr);
	free(gap->gap); free(gap->a.a); free(gap);
	return e;
}

// Using khash+ksort may be faster (identical time complexity), but unfortunately khash is 32-bit, which is not adequate
#include "kbtree.h"

typedef struct {
	int64_t x, y; // gap[x]=y
} uint128_t;

#define uint128_cmp(a, b) ((a).x - (b).x)
KBTREE_INIT(ind, uint128_t, uint128_cmp)

typedef struct {
	int64_t l0, l1, last;
	int c0, c1;
	const rld_t *e0, *e1;
	rld_t *e;
	rlditr2_t itr;
	rlditr_t itr0, itr1;
} mergeaux_t;

static void insert_tree(kbtree_t(ind) *g, uint64_t j)
{
	uint128_t z, *p;
	z.x = j;
	p = kb_getp(ind, g, &z);
	if (!p) {
		z.x = j; z.y = 1;
		kb_putp(ind, g, &z);
	} else ++p->y;
}

static kbtree_t(ind) *compute_gap_tree(const rld_t *e0, const rld_t *e1)
{
	kbtree_t(ind) *g;
	uint64_t k, l, *ok, *ol, i, j, x;
	int c = 0;
	g = kb_init(ind, KB_DEFAULT_SIZE);
	ok = alloca(8 * e0->asize);
	ol = alloca(8 * e0->asize);
	x = e1->mcnt[1];
	k = l = --x; // get the last sentinel of e1
	j = i = e0->mcnt[1] - 1; // to modify gap[j]
	insert_tree(g, j);
	for (;;) {
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
		insert_tree(g, j);
		i = j;
	}
	return g;
}

static inline void process_gap(uint128_t *p, mergeaux_t *d)
{
	//printf("gap[%lld]=%lld\n", p->x, p->y);
	dec_enc(d->e, &d->itr, d->e0, &d->itr0, &d->l0, &d->c0, p->x - d->last);
	dec_enc(d->e, &d->itr, d->e1, &d->itr1, &d->l1, &d->c1, p->y);
	d->last = p->x;
}

rld_t *fm_merge_tree(rld_t *e0, rld_t *e1, const char *fn)
{
	kbtree_t(ind) *tree;
	mergeaux_t d;
	tree = compute_gap_tree(e0, e1);
	memset(&d, 0, sizeof(mergeaux_t));
	d.last = -1;
	d.e0 = e0; d.e1 = e1;
	d.e = rld_init(e0->asize, e0->sbits, fn);
	rld_itr_init(d.e, &d.itr.itr, 0);
	d.itr.l = d.l0 = d.l1 = 0;
	d.itr.c = d.c0 = d.c1 = -1;
	rld_itr_init(e0, &d.itr0, 0);
	rld_itr_init(e1, &d.itr1, 0);
	__kb_traverse(uint128_t, tree, process_gap, &d);
	if (d.last != e0->mcnt[0] - 1)
		dec_enc(d.e, &d.itr, e0, &d.itr0, &d.l0, &d.c0, e0->mcnt[0] - 1 - d.last);
	assert(d.l0 == 0 && d.l1 == 0); // both e0 and e1 stream should be finished
	rld_enc(d.e, &d.itr.itr, d.itr.l, d.itr.c); // write the remaining symbols
	rld_enc_finish(d.e, &d.itr.itr);
	__kb_destroy(tree);
	return d.e;
}
