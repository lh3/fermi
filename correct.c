#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "utils.h"
#include "fermi.h"
#include "rld.h"
#include "kvec.h"
#include "kstring.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint32_t)

#define B_SHIFT 10
#define B_MASK ((1U<<B_SHIFT)-1)

#define MM_MAX 3
#define MM_RATIO 0.2

typedef kvec_t(uint32_t) vec32_t;

typedef struct {
	int64_t n;
	uint8_t *lock;
	vec32_t *b;
} errcorr_t;

static int g_stdout_lock;
static double g_tr, g_tc;

void seq_reverse(int l, unsigned char *s);

static inline double genpar_aux(double x, int64_t k)
{ // compute 1-(1-x)^k, where x<<1 and k is not so large
	int64_t i;
	double sum = 0., p = x, y = k;
	for (i = 1; i < k; ++i) {
		sum += p * y;
		p *= -x; y *= (double)(k - i) / (i + 1);
		if (p * y / sum < 1e-6) break;
	}
	return sum;
}

void fm_ec_genpar(int64_t n, int l, double cov, double p, int *_w, int *_T)
{
	int w, k;
	int64_t L;
	double e, qc, qe;
	L = (int64_t)((double)n * l / cov + .499);
	e = genpar_aux(p, l) * n;
	for (w = 8; w < l; ++w) {
		double q, D;
		q = genpar_aux(p, w) * (1 - p) * genpar_aux(pow(.25, w), L) * .75;
		D = genpar_aux(q, l - w) * pow(1 - p, l) * n;
		if (D < 0.0001 * e) break;
	}
	qc = (double)(l - w) / L * pow(1 - p, w + 1);
	qe = (double)(l - w) / L * (1./3.) * p * pow(1 - p, w);
	for (k = 1; k < (int)cov + 1; ++k)
		if (pow(qc, k) * pow(1 - qc, n - k) > pow(qe, k) * pow(1 - qe, n - k)) break;
	k += 2;
	fprintf(stderr, "[M::%s] HiTEC parameters for n=%ld, l=%d and c=%.1f: w_M=%d, T(w_M)=%d\n", __func__, (long)n, l, cov, w, k);
	*_w = w; *_T = k;
}

static void ec_retrieve(const rld_t *e, const fmintv_t *p, int T, kstring_t *s)
{ // get sequences descending from *p, up to bifurcation
	fmintv_t ok[6], ik = *p;
	for (;;) {
		int c, x = 0;
		fm6_extend(e, &ik, ok, 1);
		for (c = 1; c <= 4; ++c)
			if (ok[c].x[2] >= T) ++x;
		if (x == 1) {
			for (c = 1; c <= 4; ++c)
				if (ok[c].x[2] >= T) break;
			kputc(c, s);
			ik = ok[c];
		} else break;
	}
}

static void ec_save_changes(const rld_t *e, const fmintv_t *p, kstring_t *s, errcorr_t *ec, fmintv_v *stack)
{
	int c, oldl = s->l;
	fmintv_t ok[6], ik;
	size_t start = stack->n;
	ik = *p; ik.info = 0;
	kv_push(fmintv_t, *stack, ik);
	while (stack->n > start) {
		ik = kv_pop(*stack);
		s->l = oldl + (ik.info&0xffff);
		kputc(ik.info>>16, s);
		fm6_extend(e, &ik, ok, 1);
		if (ok[0].x[2]) {
			uint64_t k;
			int i, mm, l = (int)(ik.info&0xffff) + 1;
			if (l > oldl) l = oldl;
			for (i = mm = 0; i < l; ++i)
				if (s->s[i] != s->s[i + oldl]) ++mm;
			//for(i=0;i<l;++i)putchar("$ACGTN"[(int)s->s[i]]);printf(" %d %d %d\n",oldl,(int)(ik.info&0xffff)+1,mm); for (i=0;i<l;++i)putchar(s->s[i+oldl]==s->s[i]?'.':"$ACGTN"[(int)s->s[i+oldl]]);putchar('\n');
			if (mm < MM_MAX || (double)mm/l < MM_RATIO) {
				for (i = 0; i < l; ++i)
					if (s->s[i] != s->s[i + oldl]) { // an error
						uint32_t x = (uint32_t)(s->s[i] - 1)<<16 | ((ik.info&0xffff) - i);
						for (k = ok[0].x[0]; k < ok[0].x[0] + ok[0].x[2]; ++k) {
							vec32_t *b = ec->b + (k>>B_SHIFT);
							uint8_t *lock = ec->lock + (k>>B_SHIFT);
							x |= (k & B_MASK)<<18;
							while (!__sync_bool_compare_and_swap(lock, 0, 1));
							kv_push(uint32_t, *b, x);
							__sync_bool_compare_and_swap(lock, 1, 0);
						}
					}
			}
		}
		for (c = 1; c <= 5; ++c) {
			if (ok[c].x[2]) {
				ok[c].info = ((ik.info&0xffff) + 1) | c<<16;
				kv_push(fmintv_t, *stack, ok[c]);
			}
		}
	}
	s->l = oldl; s->s[s->l] = 0;
}

static void ec_collect(const rld_t *e, const fmecopt_t *opt, int len, const uint8_t *seq, errcorr_t *ec)
{
	int64_t i;
	double drop_ratio = 1. / opt->T + 1e-6;
	kstring_t str;
	fmintv_v stack;
	fmintv_t ok[6], ik;

	assert(len > 0);
	kv_init(stack);
	str.l = str.m = 0; str.s = 0;
	fm6_set_intv(e, seq[0], ik); // to get the root of the subtree to be processed
	for (i = 1; i < len; ++i) {
		fm6_extend(e, &ik, ok, 1);
		ik = ok[(int)seq[i]];
	}

	ik.info = len;
	kv_push(fmintv_t, stack, ik);
	while (stack.n) {
		int c;
		ik = kv_pop(stack);
		fm6_extend(e, &ik, ok, 1);
		if (ik.info == opt->w) { // then check and correct
			int np = 0, nn = 0; // #good/bad branches
			for (c = 1; c <= 4; ++c) {
				if (ok[c].x[2] >= opt->T) ++np;
				else if (ok[c].x[2]) ++nn;
			}
			if (np == 1 && nn) { // has one good branch and at least one bad branch(es)
				int b;
				for (b = 1; b <= 4; ++b) // base to correct to
					if (ok[b].x[2] >= opt->T) break;
				str.l = 0; kputc(b, &str); ec_retrieve(e, &ok[b], opt->T, &str); // sequence on the good branch
				for (c = 1; c <= 4; ++c) // to fix bad branch(es)
					if (ok[c].x[2] && ok[c].x[2] < opt->T && ok[c].x[2] <= opt->t && (double)ok[c].x[2] / ok[b].x[2] <= drop_ratio)
						ec_save_changes(e, &ok[c], &str, ec, &stack);
			} else if (np == 2); //	fprintf(stderr, "[E::%s] Not implemented!!!\n", __func__); // FIXME: perhaps this is necessary
		} else {
			for (c = 4; c >= 1; --c) { // FIXME: ambiguous bases are skipped
				if (ok[c].x[2] >= opt->T + 1) {
					ok[c].info = ik.info + 1;
					kv_push(fmintv_t, stack, ok[c]);
				}
			}
		}
	}

	free(stack.a); free(str.s);
}

static void ec_get_changes(const errcorr_t *ec, int64_t k, vec32_t *a)
{
	int min, max, mid, x = k&B_MASK;
	vec32_t *p = &ec->b[k>>B_SHIFT];
	if (p->n == 0) return;
	min = 0; max = p->n - 1;
	do { // binary search; linear search for small [min,max] should be faster, but the bottleneck should not be here
		mid = (min + max) / 2;
		if (x > p->a[mid]>>18) min = mid + 1;
		else max = mid - 1;
	} while (min <= max && p->a[mid]>>18 != x);
	if (p->a[mid]>>18 == x) { // we need to correct this read
		for (min = mid - 1; min >= 0 && p->a[min]>>18 == x; --min)
			kv_push(uint32_t, *a, p->a[min]);
		for (max = mid; max < p->n && p->a[max]>>18 == x; ++max)
			kv_push(uint32_t, *a, p->a[max]);
	}
}

static void ec_fix(const rld_t *e, const errcorr_t *ec, int start, int step)
{
	int64_t i, k;
	int j;
	kstring_t str, out;
	vec32_t a;

	kv_init(a);
	str.s = out.s = 0; str.l = str.m = out.l = out.m = 0;
	for (i = start<<1; i < e->mcnt[1]; i += step<<1) {
		a.n = 0;
		k = fm_retrieve(e, i + 1, &str);
		ec_get_changes(ec, k, &a);
		for (j = 0; j < a.n; ++j) // lift to the forward strand
			a.a[j] = (str.l - 1 - (a.a[j]&0xffff)) | ((3 - (a.a[j]>>16&3)) << 16);
		k = fm_retrieve(e, i, &str);
		seq_reverse(str.l, (uint8_t*)str.s);
		ec_get_changes(ec, k, &a);
		for (j = 0; j < a.n; ++j) // apply the changes
			str.s[a.a[j]&0xffff] = (a.a[j]>>16&3) + 1;
		kputc('>', &out); kputl((long)(i>>1), &out); kputc('\n', &out);
		ks_resize(&out, out.l + str.l + 2);
		for (j = 0; j < str.l; ++j)
			out.s[out.l++] = "$ACGTN"[(int)str.s[j]];
		out.s[out.l++] = '\n';
		out.s[out.l] = 0;
		if (__sync_bool_compare_and_swap(&g_stdout_lock, 0, 1)) {
			fputs(out.s, stdout);
			__sync_bool_compare_and_swap(&g_stdout_lock, 1, 0);
			out.l = 0;
		}
	}
	free(str.s); free(out.s);
}

#define MAX_DEPTH 5
#define MAX_SEQS (1<<MAX_DEPTH*2)

typedef struct {
	const rld_t *e;
	const fmecopt_t *opt;
	errcorr_t *ec;
	int n_seqs, tid;
	uint32_t seqs[MAX_SEQS];
} worker1_t;

static void *worker1(void *data)
{
	worker1_t *w = (worker1_t*)data;
	int i, j;
	for (i = 0; i < w->n_seqs; ++i) {
		uint8_t seq[MAX_DEPTH+1];
		for (j = 0; j < MAX_DEPTH; ++j)
			seq[j] = (w->seqs[i]>>(j*2)&0x3) + 1;
		ec_collect(w->e, w->opt, MAX_DEPTH, seq, w->ec);
		for (j = 0; j < MAX_DEPTH; ++j) seq[j] = "$ACGTN"[(int)seq[j]];
		seq[MAX_DEPTH] = 0;
		if (fm_verbose >= 4)
			fprintf(stderr, "[M::%s@%d] collected errors for suffix '%s' (accumulative real time: %.3f; CPU: %.3f)\n",
					__func__, w->tid, seq, realtime() - g_tr, cputime() - g_tc);
	}
	return 0;
}

typedef struct {
	const rld_t *e;
	const errcorr_t *ec;
	int start, step;
} worker2_t;

static void *worker2(void *data)
{
	worker2_t *w = (worker2_t*)data;
	ec_fix(w->e, w->ec, w->start, w->step);
	return 0;
}

int fm6_ec_correct(const rld_t *e, const fmecopt_t *opt, int n_threads)
{
	int j;
	int64_t i, avg_bucket_size;
	errcorr_t *ec;
	worker1_t *w1;
	worker2_t *w2;
	pthread_t *tid;
	pthread_attr_t attr;

	if (opt->w <= MAX_DEPTH) {
		fprintf(stderr, "[E::%s] excessively small `-w'. Please increase manually to at least %d.\n", __func__, MAX_DEPTH + 1);
		return 1;
	}
	// initialize "ec" and "tid"
	ec = calloc(1, sizeof(errcorr_t));
	ec->n = (e->mcnt[0] + B_MASK) >> B_SHIFT;
	ec->b = calloc(ec->n, sizeof(vec32_t));
	ec->lock = calloc(ec->n, 1);
	avg_bucket_size = (int64_t)((e->mcnt[0] - e->mcnt[1]) * opt->err / ec->n + .499);
	for (i = 0; i < ec->n; ++i) // preallocate memory to reduce potential locking during memory allocation
		kv_resize(uint32_t, ec->b[i], avg_bucket_size / 2);
	if (n_threads > MAX_SEQS) n_threads = MAX_SEQS;
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	// initialize and launch worker1
	g_tc = cputime(); g_tr = realtime();
	w1 = calloc(n_threads, sizeof(worker1_t));
	for (j = 0; j < n_threads; ++j)
		w1[j].ec = ec, w1[j].e = e, w1[j].opt = opt, w1[j].tid = j;
	for (i = 0, j = 0; i < MAX_SEQS; ++i) { // assign seqs
		w1[j].seqs[w1[j].n_seqs++] = i;
		if (++j == n_threads) j = 0;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker1, w1 + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(w1);
	// sort the "ec" arrays for binary search; can be multi-threaded, but should be fast enough with a single thread
	for (i = 0; i < ec->n; ++i) ks_introsort(uint32_t, ec->b[i].n, ec->b[i].a);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] collected errors in %.3f CPU seconds (%.3f wall clock)\n", __func__, cputime() - g_tc, realtime() - g_tr);

	// initialize and launch worker2
	g_tc = cputime(); g_tr = realtime();
	w2 = calloc(n_threads, sizeof(worker2_t));
	for (j = 0; j < n_threads; ++j)
		w2[j].e = e, w2[j].ec = ec, w2[j].step = n_threads, w2[j].start = j;
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker2, w2 + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(w2);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] corrected errors in %.3f CPU seconds (%.3f wall clock)\n", __func__, cputime() - g_tc, realtime() - g_tr);

	// free
	free(tid);
	for (i = 0; i < ec->n; ++i) free(ec->b[i].a);
	free(ec->b); free(ec->lock); free(ec);
	return 0;
}
