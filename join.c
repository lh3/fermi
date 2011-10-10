#include <assert.h>
#include <stdio.h>
#include <pthread.h>
#include "fermi.h"
#include "rld.h"
#include "kvec.h"
#include "kstring.h"
#include "utils.h"

static volatile int g_print_lock;
void fm_print_buffer(kstring_t *buf, volatile int *lock, int force);

static const char *g_stop_exp = "..DDDDUFBBRT";

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
			//ik.info = j - dir; kv_push(fmintv_t, *p, ik);
			ok[0].info = j - dir; kv_push(fmintv_t, *p, ok[0]);
		}
		ik = ok[c];
	}
	fm_reverse_fmivec(p); // reverse the array such that the smallest interval comes first
	return ik;
}

// requirement: s[beg..l-1] must be a full read
static int unambi_nei_for(const rld_t *e, const fmjopt_t *opt, int beg, kstring_t *s, fmintv_v *curr, fmintv_v *prev,
						  uint64_t *bits, uint64_t *term, int first, fmintv_t *lk, fmintv_t *lk0)
{
	int i, j, c, old_l = s->l, ret;
	fmintv_t ik, ok[6];
	fmintv_v *swap;
	uint64_t w[6];

	assert(s->l > opt->min_match);
	curr->n = prev->n = 0;
	// backward search for overlapping reads
	ik = overlap_intv(e, s->l - beg, (uint8_t*)s->s + beg, opt->min_match, s->l - beg - 1, 0, prev);
	if (prev->n > 0) {
		for (j = 0; j < prev->n; ++j) prev->a[j].info += beg;
		ret = prev->a[0].info; // the position with largest overlap
	} else ret = -6; // -6: no overlaps
	if (beg == 0 && first) { // the first time we call this function; test if s[0..s->l-1] contained in another read
		fm6_extend(e, &ik, ok, 1); assert(ok[0].x[2]);
		if (ik.x[2] != ok[0].x[2]) ret = -2; // the sequence is left contained
		ik = ok[0];
		fm6_extend(e, &ik, ok, 0); assert(ok[0].x[2]);
		*lk0 = ok[0]; lk0->info = (uint64_t)old_l;
		*lk = *lk0;
		if (ik.x[2] != ok[0].x[2]) ret = -3; // the sequence is right contained
		set_bits(bits, ok); // mark the read(s) has been used
	}
	if (ret < 0) return ret;
	// forward search for the forward branching test
	for (;;) {
		int c0, n_c = 0;
		memset(w, 0, 48);
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			fm6_extend(e, p, ok, 0);
			if (ok[0].x[2]) { // some reads end here
				if ((int32_t)p->info == ret && ok[0].x[2] == p->x[2]) {
					*lk = ok[0]; lk->info = (uint64_t)(old_l - ret)<<32 | (s->l - old_l); // lk keeps the terminated node
					if (bits[ok[0].x[0]>>6]>>(ok[0].x[0]&0x3f)&1) {
						ret = -10;
						set_bits(bits, ok);
						goto stop_return;
					}
					break;
				}
				set_bits(bits, ok); // mark the reads used
			}
			for (c = 1; c < 6; ++c)
				if (ok[c].x[2]) {
					if (w[c] == 0) ++n_c; // n_c keeps the number of non-zero elements in w[]
					w[c] += opt->do_dedup? 1 : ok[c].x[2];
					ok[c].info = (p->info&0xffffffffU) | (uint64_t)c<<32;
					kv_push(fmintv_t, *curr, ok[c]);
				}
		}
		if (j < prev->n || curr->n == 0) break; // found the only neighbor or cannot be extended
		if (n_c > 1) { // two or more paths
			ret = -7; goto stop_return;
			/*
			uint64_t max, sum;
			for (c0 = -1, max = sum = 0, c = 1; c < 6; ++c) {
				sum += w[c];
				if (w[c] > max) max = w[c], c0 = c;
			}
			if ((double)max / sum < 1. - opt->r || sum - max >= opt->t) { // significant conflict
				int c00 = c0;
				memset(w, 0, 48);
				for (j = 0; j < curr->n; ++j) { // count occurrences for hits longer than opt->max_match
					if (old_l - (curr->a[j].info&0xffffffff) < opt->max_match) break;
					w[curr->a[j].info>>32] += opt->do_dedup? 1 : curr->a[j].x[2];
				}
				for (c0 = -1, c = 1, sum = max = 0; c < 6; ++c) { // get the best base
					sum += w[c];
					if (w[c] > max) max = w[c], c0 = c;
				}
				if (sum == 0) c0 = c00;
				for (; j < curr->n; ++j) // continue to count the best base
					if (curr->a[j].info>>32 == c0)
						sum += opt->do_dedup? 1 : curr->a[j].x[2], max += opt->do_dedup? 1 : curr->a[j].x[2];
					else break;
				if (sum == 0 || (double)max / sum < 1. - opt->r || sum - max >= opt->t) {
					ret = -7;
					goto stop_return;
				}
			}
			for (i = j = 0; j < curr->n; ++j)
				if ((int)(curr->a[j].info>>32) == c0)
					curr->a[i++] = curr->a[j];
			ret = (int32_t)curr->a[0].info; // occasionally the nearest neighbour may be kicked out
			curr->n = i;
			*/
		} else {
			for (c0 = 1; c0 < 6; ++c0)
				if (w[c0]) break;
		}
		kputc(fm6_comp(c0), s);
		swap = curr; curr = prev; prev = swap;
	}

	// forward search for reads overlapping the extension read from the 5'-end
	overlap_intv(e, s->l, (uint8_t*)s->s, opt->max_match, ret, 1, prev);
	// backward search for backward branching test
	for (i = ret - 1; i >= beg && prev->n; --i) {
		int c00 = s->s[i], c0, n_c = 0;
		memset(w, 0, 48);
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			fm6_extend(e, p, ok, 1);
			if (ok[0].x[2]) set_bits(bits, ok);
			for (c = 1; c < 6; ++c)
				if (ok[c].x[2]) {
					if (w[c] == 0) ++n_c; // n_c keeps the number of non-zero elements in w[]
					w[c] += opt->do_dedup? 1 : ok[c].x[2];
					ok[c].info = (p->info&0xffffffffU) | (uint64_t)c<<32;
					kv_push(fmintv_t, *curr, ok[c]);
				}
		}
		if (n_c > 1) {
			ret = -8; goto stop_return;
			/*
			uint64_t max, sum;
			int i;
			for (c0 = -1, max = sum = 0, c = 1; c < 6; ++c) {
				sum += w[c];
				if (w[c] > max) max = w[c], c0 = c;
			}
			if ((double)max / sum < 1. - opt->r || sum - max >= opt->t || c0 != c00) { // ambiguous; stop extension
				ret = c0 != c00? -9 : -8;
				goto stop_return;
			}
			for (i = j = 0; j < curr->n; ++j)
				if ((int)(curr->a[j].info>>32) == c0)
					curr->a[i++] = curr->a[j];
			curr->n = i;
			*/
		}
		swap = curr; curr = prev; prev = swap;
	}
	assert(lk->x[2]);
	set_bits(bits, lk);
	if (term[lk->x[0]>>6]>>(lk->x[0]&0x3f)&1) return -11;
	return ret;

stop_return:
	s->l = old_l; s->s[s->l] = 0;
	return ret;
}

static void neighbor1(const rld_t *e, const fmjopt_t *opt, uint64_t start, uint64_t step, uint64_t *bits, uint64_t *term, fmgelem_v *elems)
{
	extern void seq_reverse(int l, unsigned char *s);
	extern void seq_revcomp6(int l, unsigned char *s);
	fmintv_v a[2];
	kstring_t s, cov;
	uint64_t x, k;

	kv_init(a[0]); kv_init(a[1]);
	s.l = s.m = cov.l = cov.m = 0; s.s = cov.s = 0;
	for (x = start<<1|1; x < e->mcnt[1]; x += step<<1) {
		int i, beg = 0, ori_len, ret[2];
		fmgelem_t z;
		fmintv_t lk[2], lk0;

		memset(lk, 0, sizeof(fmintv_t) * 2);
		k = fm_retrieve(e, x, &s);
		if (s.l <= opt->min_match) continue; // too short
		if (bits[k>>6]>>(k&0x3f)&1) continue; // the read has been used
		ks_resize(&cov, s.m);
		for (i = 0; i < s.l; ++i) cov.s[i] = '"';
		cov.l = s.l;
		{
			ori_len = s.l;
			seq_reverse(s.l, (uint8_t*)s.s);
			while ((beg = unambi_nei_for(e, opt, beg, &s, &a[0], &a[1], bits, term, 1, &lk[0], &lk0)) >= 0) {
				assert(beg);
				ks_resize(&cov, s.m);
				for (i = cov.l; i < s.l; ++i) cov.s[i] = '!';
				for (i = beg; i < s.l; ++i)
					if (cov.s[i] != '~') ++cov.s[i];
				cov.l = s.l;
			}
			lk[1] = lk0; lk[1].x[0] = lk0.x[1]; lk[1].x[1] = lk0.x[0];
			set_bit(term, lk[0].x[0]); set_bit(term, lk[0].x[1]);
			ret[0] = beg; ret[1] = 0;
		}
		if (ret[0] <= -6) { // stop due to branching or no overlaps
			beg = s.l - ori_len;
			seq_revcomp6(s.l, (uint8_t*)s.s);
			while ((beg = unambi_nei_for(e, opt, beg, &s, &a[0], &a[1], bits, term, 0, &lk[1], 0)) >= 0) {
				ks_resize(&cov, s.m);
				for (i = cov.l; i < s.l; ++i) cov.s[i] = '!';
				for (i = beg; i < s.l; ++i)
					if (cov.s[i] != '~') ++cov.s[i];
				cov.l = s.l;
			}
			set_bit(term, lk[1].x[0]); set_bit(term, lk[1].x[1]);
			ret[1] = beg;
		} else continue;
		for (i = 0; i < 2; ++i) {
			z.dir[i] = (lk[i].x[0] < lk[i].x[1])? '<' : '>';
			z.type[i] = g_stop_exp[-ret[i]];
			z.k[i] = lk[i].x[lk[i].x[0] < lk[i].x[1]? 0 : 1];
		}
		z.l = s.l;
		z.seq = malloc(z.l); z.cov = malloc(z.l);
		for (i = 0; i < s.l; ++i) {
			z.seq[i] = s.s[i] - 1;
			z.cov[i] = cov.s[i];
		}
		kv_push(fmgelem_t, *elems, z);
	}
	free(a[0].a); free(a[1].a);
	free(s.s); free(cov.s);
}

typedef struct {
	uint64_t start, step, *bits, *term;
	fmgelem_v elems;
	const rld_t *e;
	const fmjopt_t *opt;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	neighbor1(w->e, w->opt, w->start, w->step, w->bits, w->term, &w->elems);
	return 0;
}

int fm6_unambi_join(const rld_t *e, const fmjopt_t *opt, int n_threads)
{
	extern void bog_output(const fmgelem_v *elems);
	extern void bog_clean(fmgelem_v *elems);
	uint64_t *bits, *term;
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	int j;
	size_t i;
	fmgelem_v elems;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	bits = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	term = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e = e;
		ww->opt = opt;
		ww->step = n_threads;
		ww->start = j;
		ww->bits = bits;
		ww->term = term;
		kv_init(ww->elems);
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(tid); free(bits); free(term);
	elems = w[0].elems;
	for (j = 1; j < n_threads; ++j) {
		for (i = 0; i < w[j].elems.n; ++i)
			kv_push(fmgelem_t, elems, w[j].elems.a[i]);
		free(w[j].elems.a);
	}
	free(w);
	bog_clean(&elems);
	bog_output(&elems);
	free(elems.a);
	return 0;
}
