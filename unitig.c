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

static volatile int g_out_lock;

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

static fmintv_t overlap_intv(const rld_t *e, int len, const uint8_t *seq, int min, int j, int at5, fmintv_v *p, int inc_sentinel)
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
			if (inc_sentinel) {
				ok[0].info = j - dir;
				kv_push(fmintv_t, *p, ok[0]);
			} else {
				ik.info = j - dir;
				kv_push(fmintv_t, *p, ik);
			}
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
	kstring_t str;
} aux_t;

static int test_contained_right(aux_t *a, const kstring_t *s, fmintv_t *intv)
{ // read e; write prev and used
	fmintv_t ik, ok[6];
	int ret = 0;
	assert(s->l > a->min_match);
	a->a[0].n = 0;
	ik = overlap_intv(a->e, s->l, (uint8_t*)s->s, a->min_match, s->l - 1, 0, &a->a[0], 0);
	fm6_extend(a->e, &ik, ok, 1); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -1; // the sequence is left contained
	ik = ok[0];
	fm6_extend(a->e, &ik, ok, 0); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -1; // the sequence is right contained
	set_bits(a->used, ok); // mark the read(s) has been used
	*intv = ok[0];
	return ret;
}

static int try_right(aux_t *a, int beg, kstring_t *s)
{
	int ori_l = s->l, j, i, c, rbeg, is_forked = 0;
	fmintv_v *prev = &a->a[0], *curr = &a->a[1], *swap;
	fmintv_t ok[6], ok0;

	curr->n = a->nei.n = 0;
	if (prev->n == 0) { // when try_right() is called for the seed, prev is filled by test_contained_right()
		overlap_intv(a->e, s->l - beg, (uint8_t*)s->s + beg, a->min_match, s->l - beg - 1, 0, prev, 0);
		if (prev->n == 0) return -1; // no overlap
		for (j = 0; j < prev->n; ++j) prev->a[j].info += beg;
	}
	kv_resize(int, a->cat, prev->m);
	for (j = 0; j < prev->n; ++j) a->cat.a[j] = 0; // only one interval; all point to 0
	while (prev->n) {
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			if (a->cat.a[j] < 0) continue;
			fm6_extend(a->e, p, ok, 0); // forward extension
			if (ok[0].x[2] && ori_l != s->l) { // some (partial) reads end here
				fm6_extend0(a->e, &ok[0], &ok0, 1); // backward extension to look for sentinels
				if (ok0.x[2]) { // the match is bounded by sentinels - a full-length match
					if (ok[0].x[2] == p->x[2] && p->x[2] == ok0.x[2]) { // never consider a read contained in another read
						int cat = a->cat.a[j];
						assert(j == 0 || a->cat.a[j] > a->cat.a[j-1]); // otherwise not irreducible
						ok0.info = ori_l - (p->info&0xffffffffU);
						for (i = j; i < prev->n && a->cat.a[i] == cat; ++i) a->cat.a[i] = -1; // mask out other intervals of the same cat
						kv_push(fmintv_t, a->nei, ok0); // keep in the neighbor vector
						continue; // no need to go through for(c); do NOT set "used" as this neighbor may be rejected later
					} else set_bits(a->used, &ok0); // the read is contained in another read; mark it as used
				}
			} // ~if(ok[0].x[2])
			if (a->cat.a[j] < 0) continue; // no need to proceed if we have finished this path
			for (c = 1; c < 5; ++c) // collect extensible intervals
				if (ok[c].x[2]) {
					fm6_extend0(a->e, &ok[c], &ok0, 1);
					if (ok0.x[2]) { // do not extend intervals whose left end is not bounded by a sentinel
						ok[c].info = (p->info&0xfffffff0ffffffffLLU) | (uint64_t)c<<32;
						kv_push(fmintv_t, *curr, ok[c]);
					}
				}
		} // ~for(j)
		if (curr->n) { // update category
			uint32_t last, cat;
			kv_resize(int, a->cat, curr->m);
			c = curr->a[0].info>>32&0xf;
			kputc(fm6_comp(c), s);
			ks_introsort(infocmp, curr->n, curr->a);
			last = curr->a[0].info >> 32;
			a->cat.a[0] = 0;
			curr->a[0].info &= 0xffffffff;
			for (j = 1, cat = 0; j < curr->n; ++j) { // this loop recalculate cat
				if (curr->a[j].info>>32 != last) last = curr->a[j].info>>32, cat = j;
				a->cat.a[j] = cat;
				curr->a[j].info = (curr->a[j].info&0xffffffff) | (uint64_t)cat<<36;
			}
			if (cat != 0) is_forked = 1;
		}
		swap = curr; curr = prev; prev = swap; // swap curr and prev
	} // ~while(prev->n)
	if (a->nei.n == 0) return -1; // no overlap
	rbeg = ori_l - (uint32_t)a->nei.a[0].info;
	if (a->nei.n == 1 && is_forked) { // this may happen if there are contained reads
		fm6_set_intv(a->e, 0, ok0);
		for (i = rbeg; i < ori_l; ++i) {
			fm6_extend(a->e, &ok0, ok, 0);
			ok0 = ok[fm6_comp(s->s[i])];
		}
		for (i = ori_l; i < s->l; ++i) {
			int c0 = -1;
			fm6_extend(a->e, &ok0, ok, 0);
			for (c = 1, j = 0; c < 5; ++c)
				if (ok[c].x[2] && ok[c].x[0] <= a->nei.a[0].x[0] && ok[c].x[0] + ok[c].x[2] >= a->nei.a[0].x[0] + a->nei.a[0].x[2])
					++j, c0 = c;
			if (j == 0 && ok[0].x[2]) break;
			assert(j == 1);
			s->s[i] = fm6_comp(c0);
			ok0 = ok[c0];
		}
		s->l = i; s->s[s->l] = 0;
	}
	if (a->nei.n > 1) s->l = ori_l, s->s[s->l] = 0;
	return rbeg;
}

static int check_left_simple(aux_t *a, int beg, int rbeg, const kstring_t *s)
{
	fmintv_t ok[6];
	fmintv_v *prev = &a->a[0], *curr = &a->a[1], *swap;
	int i, j;

	overlap_intv(a->e, s->l, (uint8_t*)s->s, a->min_match, rbeg, 1, prev, 1);
	for (i = rbeg - 1; i >= beg; --i) {
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			fm6_extend(a->e, p, ok, 1);
			if (ok[0].x[2]) set_bits(a->used, &ok[0]); // some reads end here; they must be contained in a longer read
			if (ok[0].x[2] + ok[(int)s->s[i]].x[2] != p->x[2]) {
				return -1; // backward bifurcation
			}
			kv_push(fmintv_t, *curr, ok[(int)s->s[i]]);
		}
		swap = curr; curr = prev; prev = swap;
	} // ~for(i)
	return 0;
}

static int check_left(aux_t *a, int beg, int rbeg, const kstring_t *s)
{
	int i, ret;
	fmintv_t tmp;
	assert(a->nei.n == 1);
	ret = check_left_simple(a, beg, rbeg, s);
	if (ret == 0) return 0;
	// when ret<0, the back fork may be caused by a contained read. we have to do more to confirm this.
	tmp = a->nei.a[0]; // backup the neighbour as it will be overwritten by try_right()
	a->a[0].n = a->a[1].n = a->nei.n = 0;
	ks_resize(&a->str, s->l - rbeg + 1);
	for (i = s->l - 1, a->str.l = 0; i >= rbeg; --i)
		a->str.s[a->str.l++] = fm6_comp(s->s[i]);
	a->str.s[a->str.l] = 0;
	try_right(a, 0, &a->str);
	assert(a->nei.n >= 1);
	ret = a->nei.n > 1? -1 : 0;
	a->nei.n = 1; a->nei.a[0] = tmp; // recover the original neighbour
	return ret;
}

static void unitig_unidir(aux_t *a, kstring_t *s, kstring_t *cov, int beg0, uint64_t k0, uint64_t *end)
{ // FIXME: be careful of self-loop like a>>a or a><a
	int i, beg = beg0, rbeg, ori_l = s->l;
	while ((rbeg = try_right(a, beg, s)) >= 0) { // loop if there is at least one overlap
		uint64_t k;
		if (a->nei.n > 1) { // forward bifurcation
			set_bit(a->bend, *end);
			break;
		}
		k = a->nei.a[0].x[0];
		if (k == k0) break; // a loop like a>>b>>c>>a
		if (k == *end || a->nei.a[0].x[1] == *end) break; // a loop like a>>a or a><a
		if ((a->bend[k>>6]>>(k&0x3f)&1) || check_left(a, beg, rbeg, s) < 0) { // backward bifurcation
			set_bit(a->bend, k);
			break;
		}
		*end = a->nei.a[0].x[1];
		set_bits(a->used, &a->nei.a[0]); // successful extension
		if (cov->m < s->m) ks_resize(cov, s->m);
		cov->l = s->l; cov->s[cov->l] = 0;
		for (i = rbeg; i < ori_l; ++i) // update the coverage string
			if (cov->s[i] != '~') ++cov->s[i];
		for (i = ori_l; i < s->l; ++i) cov->s[i] = '"';
		beg = rbeg; ori_l = s->l; a->a[0].n = a->a[1].n = 0; // prepare for the next round of loop
	}
	cov->l = s->l = ori_l; s->s[ori_l] = cov->s[ori_l] = 0;
}

static int unitig1(aux_t *a, int64_t seed, kstring_t *s, kstring_t *cov, uint64_t end[2], fm128_v nei[2])
{
	extern void seq_revcomp6(int l, unsigned char *s);
	extern void seq_reverse(int l, unsigned char *s);
	fmintv_t intv0;
	int seed_len;
	fm128_t z;
	int64_t k;
	size_t i;

	nei[0].n = nei[1].n = 0;
	// retrieve the sequence pointed by seed
	k = fm_retrieve(a->e, seed, s);
	seq_reverse(s->l, (uint8_t*)s->s);
	seed_len = s->l;
	// check length, containment and if used before
	if (s->l <= a->min_match) return -1; // too short
	if (a->used[k>>6]>>(k&0x3f)&1) return -2; // used
	if (test_contained_right(a, s, &intv0) < 0) return -3; // contained; "used" is set here
	// initialize the coverage string
	if (cov->m < s->m) ks_resize(cov, s->m);
	cov->l = s->l; cov->s[cov->l] = 0;
	for (i = 0; i < cov->l; ++i) cov->s[i] = '"';
	// left-wards extension
	end[0] = intv0.x[1]; end[1] = intv0.x[0];
	if (a->a[0].n) { // no need to run this if a->a[0].n == 0
		unitig_unidir(a, s, cov, 0, intv0.x[0], &end[0]);
		for (i = 0; i < a->nei.n; ++i) {
			z.x = a->nei.a[i].x[0]; z.y = a->nei.a[i].info;
			kv_push(fm128_t, nei[0], z);
		}
	}
	// right-wards extension
	a->a[0].n = a->a[1].n = a->nei.n = 0;
	seq_revcomp6(s->l, (uint8_t*)s->s); // reverse complement for extension in the other direction
	unitig_unidir(a, s, cov, s->l - seed_len, intv0.x[1], &end[1]);
	for (i = 0; i < a->nei.n; ++i) {
		z.x = a->nei.a[i].x[0]; z.y = a->nei.a[i].info;
		kv_push(fm128_t, nei[1], z);
	}
	return 0;
}

static void write_node(kstring_t *out, long id, uint64_t kk[2], fm128_v nei[2], int l, const char *seq, const char *cov)
{
	int j, k;
	kputc('@', out); kputl(id, out);
	for (j = 0; j < 2; ++j) {
		kputc('\t', out); kputl(kk[j], out); kputc('>', out);
		for (k = 0; k != nei[j].n; ++k) {
			if (k) kputc(',', out);
			kputl(nei[j].a[k].x, out); kputc(':', out); kputw(nei[j].a[k].y, out);
		}
		if (nei[j].n == 0) kputc('.', out);
	}
	kputc('\n', out);
	ks_resize(out, out->l + 2 * l + 5);
	for (j = 0; j < l; ++j) out->s[out->l++] = "ACGT"[(int)seq[j] - 1];
	kputsn("\n+\n", 3, out);
	kputsn(cov, l, out);
	kputc('\n', out);
}

static void unitig_core(const rld_t *e, int min_match, int64_t start, int64_t end, uint64_t *used, uint64_t *bend, uint64_t *visited)
{
	extern void fm_print_buffer(kstring_t *buf, volatile int *lock, int force);
	uint64_t i, k[2];
	aux_t a;
	kstring_t str, cov, out;
	fm128_v nei[2];

	assert((start&1) == 0 && (end&1) == 0);
	// initialize aux_t and all the vectors
	a.str.l = a.str.m = str.l = str.m = cov.l = cov.m = out.l = out.m = 0; str.s = cov.s = out.s = 0;
	a.e = e; a.min_match = min_match; a.used = used; a.bend = bend;
	kv_init(a.a[0]); kv_init(a.a[1]); kv_init(a.nei); kv_init(a.cat);
	kv_init(nei[0]); kv_init(nei[1]);
	// the core loop
	for (i = start|1; i < end; i += 2) {
		if (unitig1(&a, i, &str, &cov, k, nei) >= 0) { // then we keep the unitig
			uint64_t *p[2], x[2];
			p[0] = visited + (k[0]>>6); x[0] = 1LLU<<(k[0]&0x3f);
			p[1] = visited + (k[1]>>6); x[1] = 1LLU<<(k[1]&0x3f);
			if ((__sync_fetch_and_or(p[0], x[0])&x[0]) || (__sync_fetch_and_or(p[1], x[1])&x[1])) continue; // NOT always working
			write_node(&out, i, k, nei, str.l, str.s, cov.s);
			fm_print_buffer(&out, &g_out_lock, 0);
		}
	}
	fm_print_buffer(&out, &g_out_lock, 1);
	free(a.a[0].a); free(a.a[1].a); free(a.nei.a); free(a.cat.a);
	free(nei[0].a); free(nei[1].a);
	free(a.str.s); free(str.s); free(cov.s); free(out.s);
}

typedef struct {
	uint64_t start, end, *used, *bend, *visited;
	const rld_t *e;
	int min_match;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	unitig_core(w->e, w->min_match, w->start, w->end, w->used, w->bend, w->visited);
	return 0;
}

int fm6_unitig(const rld_t *e, int min_match, int n_threads)
{
	uint64_t *used, *bend, *visited, rest = e->mcnt[1];
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	int j;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	used    = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	bend    = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	visited = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	assert(e->mcnt[1] >= n_threads * 2);
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e = e;
		ww->min_match = min_match;
		ww->start = (e->mcnt[1] - rest) / 2 * 2;
		ww->end = ww->start + rest / (n_threads - j) / 2 * 2;
		rest -= ww->end - ww->start;
		ww->used = used; ww->bend = bend; ww->visited = visited;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(tid); free(used); free(bend); free(visited); free(w);
	return 0;
}
