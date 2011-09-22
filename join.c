#include <assert.h>
#include <stdio.h>
#include <pthread.h>
#include "fermi.h"
#include "rld.h"
#include "kvec.h"
#include "kstring.h"
#include "utils.h"

static uint64_t g_cnt;
static int g_print_lock;

static inline void set_bit(uint64_t *bits, uint64_t x)
{
	uint64_t *p = bits + (x>>6);
	uint64_t z = 1LLU<<(x&0x3f), y;
	y = __sync_fetch_and_or(p, z);
	if ((y & z) == 0) {
		__sync_add_and_fetch(&g_cnt, 1);
	}
}

static inline void set_bits(uint64_t *bits, const fmintv_t *p)
{
	uint64_t k;
	for (k = 0; k < p->x[2]; ++k) {
		set_bit(bits, p->x[0] + k);
		set_bit(bits, p->x[1] + k);
	}
}

// requirement: s[beg..l-1] must be a full read
static int unambi_nei_for(const rld_t *e, int min, int beg, kstring_t *s, fmintv_v *curr, fmintv_v *prev, uint64_t *bits)
{
	extern fmintv_t fm6_overlap_intv(const rld_t *e, int len, const uint8_t *seq, int min, int j, int at5, fmintv_v *p);
	int i, j, c, old_l = s->l, ret;
	fmintv_t ik, ok[6];
	fmintv_v *swap;

	curr->n = prev->n = 0;
	// backward search for overlapping reads
	ik = fm6_overlap_intv(e, s->l - beg, (uint8_t*)s->s + beg, min, s->l - beg - 1, 0, prev);
	assert(prev->n || (int)s->l < min);
	if (prev->n > 0) {
		for (j = 0; j < prev->n; ++j) prev->a[j].info += beg;
		ret = prev->a[0].info; // the position with largest overlap
	} else ret = -1; // read is too short
	// test if s[beg..s->l-1] contained in another read
	fm6_extend(e, &ik, ok, 1); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -2; // the sequence is left contained
	ik = ok[0];
	fm6_extend(e, &ik, ok, 0); assert(ok[0].x[2]);
	if (ik.x[2] != ok[0].x[2]) ret = -3; // the sequence is right contained
	set_bits(bits, ok); // mark the read(s) has been used
	if (ret < 0) return ret;
	// forward search for the forward branching test
	while (prev->n) {
		int c0 = -1;
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			fm6_extend(e, p, ok, 0);
			if (ok[0].x[2]) {
				set_bits(bits, ok);
				if ((int)p->info == ret && ok[0].x[2] == p->x[2]) break;
			}
			if (c0 == -1) {
				for (c = 1; c < 6; ++c) if (ok[c].x[2]) break;
				if (c == 6) continue;
				c0 = c;
			}
			if (ok[c0].x[2] + ok[0].x[2] < p->x[2]) return -4;
			if (ok[c0].x[2] && (curr->n == 0 || ok[c0].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c0].info = p->info;
				kv_push(fmintv_t, *curr, ok[c0]);
			}
		}
		if (j < prev->n) break;
		kputc(fm6_comp(c0), s);
		swap = curr; curr = prev; prev = swap;
	}
	for (i = 0; i < s->l; ++i) putchar("$ACGTN"[(int)s->s[i]]); putchar('\n');

	// forward search for reads overlapping the extension read from the 5'-end
	fm6_overlap_intv(e, s->l, (uint8_t*)s->s, min, ret, 1, prev);
	printf("ret=%d, len=%d, prev->n=%d\n", (int)ret, (int)s->l, (int)prev->n);
	// backward search for backward branching test
	for (i = ret - 1; i >= 0 && prev->n; --i) {
		c = s->s[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fm6_extend(e, &prev->a[j], ok, 1);
			if (ok[0].x[2]) set_bits(bits, ok);
			if (ok[c].x[2] + ok[0].x[2] != prev->a[j].x[2]) { // branching
				s->l = old_l;
				return -5; // backward branching
			}
			if (ok[c].x[2] && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]))
				kv_push(fmintv_t, *curr, ok[c]);
		}
		swap = curr; curr = prev; prev = swap;
	}
	printf("final: ret=%d, len=%d\n", (int)ret, (int)s->l);
	return ret;
}

static void neighbor1(const rld_t *e, int min, uint64_t start, uint64_t step, uint64_t *bits, FILE *fp)
{
	extern void seq_reverse(int l, unsigned char *s);
	extern void seq_revcomp6(int l, unsigned char *s);
	fmintv_v a[2];
	kstring_t s, out;
	uint64_t x, k;

	kv_init(a[0]); kv_init(a[1]);
	s.l = s.m = 0; s.s = 0;
	out.l = out.m = 0; out.s = 0;
	for (x = start; x < e->mcnt[1]; x += step) {
		int i, beg = 0, ori_len;
		k = fm_retrieve(e, x, &s);
		if (bits[k>>6]>>(k&0x3f)&1) continue; // the read has been used
		ori_len = s.l;
		seq_reverse(s.l, (uint8_t*)s.s);
		while ((beg = unambi_nei_for(e, min, beg, &s, &a[0], &a[1], bits)) >= 0);
		if (beg <= -4) { // stop due to branching
			beg = s.l - ori_len;
			seq_revcomp6(s.l, (uint8_t*)s.s);
			while ((beg = unambi_nei_for(e, min, beg, &s, &a[0], &a[1], bits)) >= 0);
		}
		kputc('>', &out); kputl((long)x, &out); kputc('\n', &out);
		for (i = 0; i < s.l; ++i)
			kputc("$ACGTN"[(int)s.s[i]], &out);
		kputc('\n', &out);
		if (__sync_bool_compare_and_swap(&g_print_lock, 0, 1)) {
			fputs(out.s, fp);
			out.l = 0;
			__sync_bool_compare_and_swap(&g_print_lock, 1, 0);
		}
	}
	while (__sync_bool_compare_and_swap(&g_print_lock, 0, 1) == 0); // busy waiting
	fputs(out.s, fp);
	__sync_bool_compare_and_swap(&g_print_lock, 1, 0);
	free(a[0].a); free(a[1].a);
	free(s.s); free(out.s);
}

typedef struct {
	int min;
	uint64_t start, step, *bits;
	const rld_t *e;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	neighbor1(w->e, w->min, w->start, w->step, w->bits, stdout);
	return 0;
}

int fm6_unambi_join(const rld_t *e, int min, int n_threads)
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
		ww->step = n_threads;
		ww->start = j;
		ww->bits = bits;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(w); free(tid); free(bits);
	return 0;
}
