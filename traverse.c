#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "fermi.h"
#include "rld.h"
#include "kvec.h"

#define B_SHIFT 14
#define B_MASK ((1U<<B_SHIFT)-1)

typedef struct {
	int64_t n, m;
	uint32_t *a;
} ec_bucket_t;

typedef struct {
	int64_t n;
	ec_bucket_t *b;
} errcorr_t;

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

errcorr_t *fm6_ec_collect(const rld_t *e, const fmecopt_t *opt, int len, const uint8_t *seq)
{
	errcorr_t *ec;
	int64_t i, tot = 0, *cnt;
	double drop_ratio = 1. / (opt->T + 1);
	fmintv_v stack;
	fmintv_t ok[6], ik;

	assert(len > 0);
	kv_init(stack);
	ec = calloc(1, sizeof(errcorr_t));
	ec->n = (e->mcnt[0] + B_MASK) >> B_SHIFT;
	ec->b = calloc(ec->n, sizeof(ec_bucket_t));
	cnt = alloca(8 * (opt->depth + 1));
	memset(cnt, 0, 8 * (opt->depth + 1));

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
		++cnt[ik.info];
		++tot;
		if (ik.info == opt->depth) {
		} else {
			for (c = 4; c >= 1; --c) { // FIXME: ambiguous bases
				if (ok[c].x[2] >= opt->T + 1) {
					ok[c].info = ik.info + 1;
					kv_push(fmintv_t, stack, ok[c]);
				}
			}
		}
	}

//	for (i = 1; i <= opt->depth; ++i) fprintf(stderr, "%lld, %lld, %g\n", i, cnt[i], (double)cnt[i] / pow(4, i - len));
	fprintf(stderr, "%lld\n", tot);
	free(stack.a);
	return ec;
}

int fm6_ec_correct(const rld_t *e, const fmecopt_t *opt, int n_threads)
{
	uint8_t c;
	fm_ec_genpar(e->mcnt[1]/2, (int)((double)e->mcnt[0] / e->mcnt[1] + 0.5), 30.0, 0.01);
	for (c = 1; c <= 4; ++c)
		fm6_ec_collect(e, opt, 1, &c);
	return 0;
}
