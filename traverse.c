#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "fermi.h"
#include "rld.h"
#include "kvec.h"
#include "kstring.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint32_t)

#define B_SHIFT 10
#define B_MASK ((1U<<B_SHIFT)-1)

typedef kvec_t(uint32_t) vec32_t;

typedef struct {
	int64_t n;
	vec32_t *b;
} errcorr_t;

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

void fm_ec_genpar(int64_t n, int l, double cov, double p)
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
}

static void save_fix(const rld_t *e, const fmintv_t *p, int b, errcorr_t *ec, fmintv_v *stack)
{
	int c;
	fmintv_t ok[6], ik;
	size_t start = stack->n;
	ik = *p; ik.info = 0;
	kv_push(fmintv_t, *stack, ik);
	while (stack->n > start) {
		ik = kv_pop(*stack);
		fm6_extend(e, &ik, ok, 1);
		if (ok[0].x[2]) {
			uint64_t k;
			for (k = ok[0].x[0]; k < ok[0].x[0] + ok[0].x[2]; ++k) {
				uint32_t x = (uint32_t)b<<16 | (k & B_MASK)<<18 | ik.info;
				vec32_t *b = ec->b + (k>>B_SHIFT);
				kv_push(uint32_t, *b, x);
			}
		}
		for (c = 1; c <= 5; ++c) {
			if (ok[c].x[2]) {
				ok[c].info = ik.info + 1;
				kv_push(fmintv_t, *stack, ok[c]);
			}
		}
	}
}

static void ec_collect(const rld_t *e, const fmecopt_t *opt, int len, const uint8_t *seq, errcorr_t *ec)
{
	int64_t i;
	double drop_ratio = 1. / opt->T + 1e-6;
	fmintv_v stack;
	fmintv_t ok[6], ik;

	assert(len > 0);
	kv_init(stack);
	fm6_set_intv(e, seq[0], ik);
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
		if (ik.info == opt->depth) {
			int np = 0, nn = 0;
			for (c = 1; c <= 4; ++c) {
				if (ok[c].x[2] >= opt->T) ++np;
				else if (ok[c].x[2]) ++nn;
			}
			if (np == 1) {
				int b;
				for (b = 1; b <= 4; ++b) // base to correct to
					if (ok[b].x[2] >= opt->T) break;
				for (c = 1; c <= 4; ++c)
					if (ok[c].x[2] && ok[c].x[2] < opt->T && ok[c].x[2] <= opt->t && (double)ok[c].x[2] / ok[b].x[2] <= drop_ratio)
						save_fix(e, &ok[c], b - 1, ec, &stack);
			} else if (np == 2) { // FIXME: not implemented
			}
		} else {
			for (c = 4; c >= 1; --c) { // FIXME: ambiguous bases
				if (ok[c].x[2] >= opt->T + 1) {
					ok[c].info = ik.info + 1;
					kv_push(fmintv_t, stack, ok[c]);
				}
			}
		}
	}

	free(stack.a);
}

static void ec_get_changes(const errcorr_t *ec, int64_t k, vec32_t *a)
{
	int min, max, mid, x = k&B_MASK;
	vec32_t *p = &ec->b[k>>B_SHIFT];
	if (p->n == 0) return;
	min = 0; max = p->n - 1;
	do { // binary search
		mid = (min + max) / 2;
		if (x > p->a[mid]>>18) min = mid + 1;
		else max = mid - 1;
	} while (min < max && p->a[mid]>>18 != x);
	if (p->a[mid]>>18 == x) {
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
	start = start >> 1 << 1;
	for (i = start; i < e->mcnt[1]; i += step<<1) {
		a.n = 0;
		k = fm_retrieve(e, i + 1, &str);
		ec_get_changes(ec, k, &a);
		for (j = 0; j < a.n; ++j) // lift to the forward strand
			a.a[j] = (str.l - 1 - (a.a[j]&0xffff)) | ((3 - (a.a[j]>>16&3)) << 16);
		k = fm_retrieve(e, i, &str);
		seq_reverse(str.l, (uint8_t*)str.s);
		ec_get_changes(ec, k, &a);
		for (j = 0; j < a.n; ++j)
			str.s[a.a[j]&0xffff] = (a.a[j]>>16&3) + 1;
		kputc('>', &out); kputl((long)i, &out); kputc('\n', &out);
		ks_resize(&out, out.l + str.l + 2);
		for (j = 0; j < str.l; ++j)
			out.s[out.l++] = "$ACGTN"[(int)str.s[j]];
		out.s[out.l++] = '\n';
		out.s[out.l] = 0;
		{
			fputs(out.s, stdout);
			out.l = 0;
		}
	}
	free(str.s); free(out.s);
}

int fm6_ec_correct(const rld_t *e, const fmecopt_t *opt, int n_threads)
{
	uint8_t c;
	int64_t i;
	errcorr_t *ec;
	ec = calloc(1, sizeof(errcorr_t));
	ec->n = (e->mcnt[0] + B_MASK) >> B_SHIFT;
	ec->b = calloc(ec->n, sizeof(vec32_t));
	fm_ec_genpar(e->mcnt[1]/2, (int)((double)e->mcnt[0] / e->mcnt[1] + 0.5), 30.0, 0.01);
	for (c = 1; c <= 4; ++c)
		ec_collect(e, opt, 1, &c, ec);
	for (i = 0; i < ec->n; ++i)
		ks_introsort(uint32_t, ec->b[i].n, ec->b[i].a);
	ec_fix(e, ec, 0, 1);
	return 0;
}
