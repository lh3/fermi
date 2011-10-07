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
static int unambi_nei_for(const rld_t *e, const fmjopt_t *opt, int beg, kstring_t *s, fmintv_v *curr, fmintv_v *prev, uint64_t *bits, int first)
{
	int i, j, c, old_l = s->l, ret;
	fmintv_t ik, ok[6];
	fmintv_v *swap;
	uint64_t w[6];

	curr->n = prev->n = 0;
	// backward search for overlapping reads
	ik = overlap_intv(e, s->l - beg, (uint8_t*)s->s + beg, opt->min_match, s->l - beg - 1, 0, prev);
	//for (i = 0, c = 0; i < prev->n; ++i) c += prev->a[i].x[2]; fprintf(stderr, "Total: %d\n", c);
	if (prev->n > 0) {
		for (j = 0; j < prev->n; ++j) prev->a[j].info += beg;
		ret = prev->a[0].info; // the position with largest overlap
	} else ret = s->l - beg <= opt->min_match? -1 : -6; // -1: too short; -6: no overlaps
	if (beg == 0 && first) { // test if s[beg..s->l-1] contained in another read
		fm6_extend(e, &ik, ok, 1); assert(ok[0].x[2]);
		if (ik.x[2] != ok[0].x[2]) ret = -2; // the sequence is left contained
		ik = ok[0];
		fm6_extend(e, &ik, ok, 0); assert(ok[0].x[2]);
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
					if (bits[ok[0].x[0]>>6]>>(ok[0].x[0]&0x3f)&1) ret = -10;
					else set_bits(bits, ok);
					if (ret < 0) goto stop_return;
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
		if (curr->n == 0) break; // cannot be extended
		if (j < prev->n) break; // found the only neighbor
		if (n_c > 1) { // two or more paths
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
		} else {
			for (c0 = 1; c0 < 6; ++c0)
				if (w[c0]) break;
		}
		kputc(fm6_comp(c0), s);
		swap = curr; curr = prev; prev = swap;
	}
	//for (i = 0; i < s->l; ++i) putchar("$ACGTN"[(int)s->s[i]]); putchar('\n');

	// forward search for reads overlapping the extension read from the 5'-end
	overlap_intv(e, s->l, (uint8_t*)s->s, opt->max_match, ret, 1, prev);
	//printf("ret=%d, len=%d, prev->n=%d, %d\n", (int)ret, (int)s->l, (int)prev->n, min);
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
		}
		swap = curr; curr = prev; prev = swap;
	}
	//printf("final: ret=%d, len=%d\n", (int)ret, (int)s->l);
	return ret;

stop_return:
	s->l = old_l; s->s[s->l] = 0;
	return ret;
}

static void neighbor1(const rld_t *e, const fmjopt_t *opt, uint64_t start, uint64_t step, uint64_t *bits, FILE *fp)
{
	extern void seq_reverse(int l, unsigned char *s);
	extern void seq_revcomp6(int l, unsigned char *s);
	fmintv_v a[2];
	kstring_t s, out;
	uint64_t x, k;

	kv_init(a[0]); kv_init(a[1]);
	s.l = s.m = 0; s.s = 0;
	out.l = out.m = 0; out.s = 0;
	for (x = start<<1|1; x < e->mcnt[1]; x += step<<1) {
		int i, beg = 0, ori_len, ret1, left_cnt = 0, rght_cnt = 0;
		k = fm_retrieve(e, x, &s);
		if (bits[k>>6]>>(k&0x3f)&1) continue; // the read has been used
		ori_len = s.l;
		seq_reverse(s.l, (uint8_t*)s.s);
		while ((beg = unambi_nei_for(e, opt, beg, &s, &a[0], &a[1], bits, 1)) >= 0) ++left_cnt;
		if ((ret1 = beg) <= -6) { // stop due to branching or no overlaps
			beg = s.l - ori_len;
			seq_revcomp6(s.l, (uint8_t*)s.s);
			while ((beg = unambi_nei_for(e, opt, beg, &s, &a[0], &a[1], bits, 0)) >= 0) ++rght_cnt;
		} else continue;
		if (left_cnt == 0 && rght_cnt == 0) continue;
		kputc('>', &out); kputl((long)x, &out); kputc(' ', &out); kputw(ret1, &out); kputw(beg, &out); kputc(' ', &out);
		kputw(left_cnt, &out); kputc(' ', &out); kputw(rght_cnt, &out); kputc('\n', &out);
		for (i = 0; i < s.l; ++i)
			kputc("$ACGTN"[(int)s.s[i]], &out);
		kputc('\n', &out);
		fm_print_buffer(&out, &g_print_lock, 0);
	}
	fm_print_buffer(&out, &g_print_lock, 1);
	free(a[0].a); free(a[1].a);
	free(s.s); free(out.s);
}

typedef struct {
	uint64_t start, step, *bits;
	const rld_t *e;
	const fmjopt_t *opt;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	neighbor1(w->e, w->opt, w->start, w->step, w->bits, stdout);
	return 0;
}

int fm6_unambi_join(const rld_t *e, const fmjopt_t *opt, int n_threads)
{
	uint64_t *bits;
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	int j;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	bits = (uint64_t*)xcalloc((e->mcnt[1] + 63)/64, 8);
	for (j = 0; j < n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e = e;
		ww->opt = opt;
		ww->step = n_threads;
		ww->start = j;
		ww->bits = bits;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(w); free(tid); free(bits);
	return 0;
}
