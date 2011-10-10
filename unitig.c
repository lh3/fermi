#include <assert.h>
#include <pthread.h>
#include <string.h>
#include "fermi.h"
#include "rld.h"
#include "kvec.h"
#include "kstring.h"
#include "utils.h"

#define info_lt(a, b) ((a).info < (b).info)

#include "ksort.h"
KSORT_INIT(infocmp, fmintv_t, info_lt)

static inline void set_bit(uint64_t *bits, uint64_t x)
{
	uint64_t *p = bits + (x>>6);
	uint64_t z = 1LLU<<(x&0x3f);
	__sync_fetch_and_or(p, z);
}

static inline void set_bits(uint64_t *bits, const fmintv_t *p)
{
	uint64_t k;
	for (k = 0; k < p->x[2]; ++k) {
		set_bit(bits, p->x[0] + k);
		set_bit(bits, p->x[1] + k);
	}
}

static fmintv_t overlap_intv(const rld_t *e, int len, const uint8_t *seq, int min, int j, int at5, fmintv_v *p)
{ // requirement: seq[j] matches the end of a read
	extern void fm_reverse_fmivec(fmintv_v *p);
	int c, depth, dir, end;
	fmintv_t ik, ok[6];
	p->n = 0;
	dir = at5? 1 : -1; // at5 is true iff we start from the 5'-end of a read
	end = at5? len : -1;
	c = seq[j];
	fm6_set_intv(e, c, ik);
	for (depth = 1, j += dir; j != end; j += dir, ++depth) {
		c = at5? fm6_comp(seq[j]) : seq[j];
		fm6_extend(e, &ik, ok, !at5);
		if (!ok[c].x[2]) break; // cannot be extended
		if (depth >= min && ok[0].x[2]) {
			ok[0].info = j - dir;
			kv_push(fmintv_t, *p, ok[0]);
		}
		ik = ok[c];
	}
	fm_reverse_fmivec(p); // reverse the array such that the smallest interval comes first
	return ik;
}

typedef struct {
	const rld_t *e;
	int min_match;
	fmintv_v a[2], nei;
	kvec_t(int) cat;
	uint64_t *used, *bend;
} aux_t;

static int test_contained_right(aux_t *a, const kstring_t *s, fmintv_t *intv)
{ // read e; write prev and used
	fmintv_t ik, ok[6];
	int ret = 0;
	assert(s->l > a->min_match);
	a->a[0].n = 0;
	ik = overlap_intv(a->e, s->l, (uint8_t*)s->s, a->min_match, s->l - 1, 0, &a->a[0]);
	fm6_extend(a->e, &ik, ok, 1); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -1; // the sequence is left contained
	ik = ok[0];
	fm6_extend(a->e, &ik, ok, 0); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -1; // the sequence is right contained
	set_bits(a->used, ok); // mark the read(s) has been used
	*intv = ok[0];
	return ret;
}

static int extend_right(aux_t *a, int beg, kstring_t *s)
{
	int ori_l = s->l, j, i, c, rbeg;
	fmintv_v *prev = &a->a[0], *curr = &a->a[1], *swap;
	fmintv_t ok[6];

	curr->n = a->nei.n = 0;
	if (prev->n == 0) { // when extend_right() is called for the seed, prev is filled by test_contained_right()
		overlap_intv(a->e, s->l - beg, (uint8_t*)s->s + beg, a->min_match, s->l - beg - 1, 0, prev);
		if (prev->n == 0) return -1; // no overlap
		for (j = 0; j < prev->n; ++j) prev->a[j].info += beg;
	}
	kv_resize(int, a->cat, prev->m);
	for (j = 0; j < prev->n; ++j) a->cat.a[j] = 0; // only one interval; all point to 0
	rbeg = prev->a[0].info&0xffffffffU; // read start
	for (;;) {
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			if (a->cat.a[j] < 0) continue;
			fm6_extend(a->e, p, ok, 0);
			if (ok[0].x[2]) { // some reads end here
				if ((int32_t)p->info == (int32_t)prev->a[a->cat.a[j]].info && ok[0].x[2] == p->x[2]) { // found the neighbor
					ok[0].info = s->l - ori_l;
					kv_push(fmintv_t, a->nei, ok[0]); // keep in the neighbor vector
					// mask out other intervals of the same cat(egory)
					for (i = j + 1; i < prev->n && a->cat.a[i] == a->cat.a[j]; ++i) a->cat.a[i] = -1;
					a->cat.a[j] = -1;
					continue; // no need to go through for(c); do NOT set "used" as this neighbor may be rejected later
				}
				set_bits(a->used, ok); // the read is contained in another read; mark it as used
			}
			for (c = 1; c < 5; ++c) // collect extensible intervals
				if (ok[c].x[2]) {
					ok[c].info = (p->info&0xfffffff0ffffffffLLU) | (uint64_t)c<<32;
					kv_push(fmintv_t, *curr, ok[c]);
				}
		}
		if (curr->n == 0) break;
		{ // update category
			uint32_t last, cat;
			c = curr->a[0].info>>32&0xf;
			kputc(fm6_comp(c), s);
			ks_introsort(infocmp, curr->n, curr->a);
			last = curr->a[0].info >> 32;
			a->cat.a[0] = 0;
			for (j = 1, cat = 0; j < prev->n; ++j) {
				if (curr->a[j].info>>32 != last)
					last = curr->a[j].info>>32, cat = j, curr->a[j].info = (curr->a[j].info&0xffffffff) | (uint64_t)cat<<32;
				a->cat.a[j] = cat;
			}
		}
		swap = curr; curr = prev; prev = swap;
	}
	if (a->nei.n > 1) s->l = ori_l, s->s[s->l] = 0;
	return rbeg;
}

static int check_left(aux_t *a, int beg, int rbeg, const kstring_t *s)
{
	fmintv_t ok[6];
	fmintv_v *prev = &a->a[0], *curr = &a->a[1], *swap;
	int i, j;

	overlap_intv(a->e, s->l, (uint8_t*)s->s, a->min_match, rbeg, 1, prev);
	for (i = rbeg - 1; i >= beg; --i) {
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			fm6_extend(a->e, p, ok, 1);
			if (ok[0].x[2]) set_bits(a->used, ok); // some reads end here; they must be contained in a longer read
			if (ok[0].x[2] + ok[(int)s->s[i]].x[2] != p->x[2]) return -1; // backward bifurcation
			kv_push(fmintv_t, *curr, ok[(int)s->s[i]]);
		}
		swap = curr; curr = prev; prev = swap;
	}
	return 0;
}

static void unitig_unidir(aux_t *a, kstring_t *s, uint64_t k0, uint64_t *end)
{ // FIXME: be careful of self-loop like a>>a or a><a
	int beg = 0, rbeg, old_l = s->l;
	while ((rbeg = extend_right(a, beg, s)) >= 0) { // loop if there is at least one overlap
		uint64_t k = a->nei.a[0].x[0];
		if (a->nei.n > 1) break; // forward bifurcation
		if (k == k0) break; // a loop like a>>b>>c>>a
		if (k == *end || a->nei.a[0].x[1] == *end) break; // a loop like a>>a or a><a
		if ((a->bend[k>>6]>>(k&0x3f)&1) || check_left(a, beg, rbeg, s) < 0) { // backward bifurcation
			s->l = old_l; s->s[old_l] = 0; // revert to the sequence before extension
			set_bit(a->bend, k);
			return;
		}
		*end = a->nei.a[0].x[1];
		set_bits(a->used, &a->nei.a[0]); // successful extension
		beg = rbeg; old_l = s->l; a->a[0].n = 0; // prepare for the next round of loop
	}
}

static int unitig1(aux_t *a, int64_t seed, kstring_t *s, uint64_t end[2], fm64_v nei[2])
{
	extern void seq_revcomp6(int l, unsigned char *s);
	extern void seq_reverse(int l, unsigned char *s);
	int64_t k;
	size_t i;
	fmintv_t intv0;

	nei[0].n = nei[1].n = 0;
	k = fm_retrieve(a->e, seed, s);
	seq_reverse(s->l, (uint8_t*)s->s);
	if (s->l <= a->min_match) return -1; // too short
	if (a->used[k>>6]>>(k&0x3f)&1) return -2; // used
	if (test_contained_right(a, s, &intv0) < 0) return -3; // contained; "used" is set here
	end[0] = intv0.x[1]; end[1] = intv0.x[0];
	unitig_unidir(a, s, intv0.x[0], &end[0]);
	for (i = 0; i < a->nei.n; ++i) kv_push(uint64_t, nei[0], a->nei.a[i].x[0]);
	seq_revcomp6(s->l, (uint8_t*)s->s);
	unitig_unidir(a, s, intv0.x[1], &end[1]);
	for (i = 0; i < a->nei.n; ++i) kv_push(uint64_t, nei[1], a->nei.a[i].x[0]);
	return 0;
}

static void unitig_core(const rld_t *e, int min_match, int64_t start, int64_t step, uint64_t *used, uint64_t *bend, fmnode_v *nodes)
{
	uint64_t i, end[2];
	int j;
	aux_t a;
	kstring_t str;
	fm64_v nei[2];

	str.l = str.m = 0; str.s = 0;
	a.e = e; a.min_match = min_match; a.used = used; a.bend = bend;
	kv_init(a.a[0]); kv_init(a.a[1]); kv_init(a.nei); kv_init(a.cat);
	kv_init(nei[0]); kv_init(nei[1]);
	for (i = start<<1|1; i < e->mcnt[1]; i += step<<1) {
		for (j = 0; j < 4; j += step) {
			if (unitig1(&a, i + j, &str, end, nei) >= 0) {
				fmnode_t *p;
				kv_pushp(fmnode_t, *nodes, &p);
				p->k[0] = end[0]; kv_init(p->nei[0]); kv_copy(uint64_t, p->nei[0], nei[0]);
				p->k[1] = end[1]; kv_init(p->nei[1]); kv_copy(uint64_t, p->nei[1], nei[1]);
				p->l = str.l;
				p->seq = calloc(p->l, 1);
				memcpy(p->seq, str.s, p->l);
			}
		}
	}
	free(a.a[0].a); free(a.a[1].a); free(a.nei.a); free(a.cat.a);
	free(nei[0].a); free(nei[1].a);
}

typedef struct {
	uint64_t start, step, *used, *bend;
	const rld_t *e;
	int min_match;
	fmnode_v nodes;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	unitig_core(w->e, w->min_match, w->start, w->step, w->used, w->bend, &w->nodes);
	return 0;
}

int fm6_unitig(const rld_t *e, int min_match, int n_threads)
{
	extern void msg_print(const fmnode_v *nodes);
	uint64_t *used, *bend;
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	fmnode_v nodes;
	int j;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	used  = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	bend = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e = e;
		ww->min_match = min_match;
		ww->step = n_threads * 2;
		ww->start = j;
		ww->used = used;
		ww->bend = bend;
		kv_init(ww->nodes);
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(tid); free(used); free(bend);
	kv_init(nodes);
	if (n_threads > 1) {
		size_t i;
		for (j = 0; j < n_threads; ++j) {
			for (i = 0; i < w[j].nodes.n; ++i)
				kv_push(fmnode_t, nodes, w[j].nodes.a[i]);
			free(w[j].nodes.a);
		}
	} else nodes = w[0].nodes;
	msg_print(&nodes);
	free(w);
	return 0;
}
