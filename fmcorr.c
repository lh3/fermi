#include "fermi.h"
#include "rld.h"
#include "kvec.h"

typedef struct {
	int min_l, min_occ;
	double ratio_factor;
	double diff_factor;
	int for_qthres;
} fmcorr_opt_t;

typedef struct {
	fmintv_t k[6];
} fmintv6_t;

typedef struct { size_t n, m; fmintv6_t *a; } fmintv6_v;

static inline int ec2_score(const fmcorr_opt_t *opt, const fmintv_t k[6], int c, int *ec)
{
	int i;
	uint64_t sum = 0;
	for (i = 1; i <= 4; ++i) sum += k[i];
}

int fm_ec2_core(const fmcorr_opt_t *opt, const rld_t *e, int l_seq, uint8_t *seq, uint8_t *qual, int x, fmintv_v *prev, fmintv_v *curr, fmintv6_v *tmp)
{ // this function is similar to fm6_smem1_core()
	int i, j, c, ret;
	fmintv_t ik;
	fmintv_v *swap;
	fmintv6_t ok;
	
	fm6_set_intv(e, seq[x], ik);
	ik.info = x + 1;

	for (i = x + 1; i < l_seq; ++i) { // forward search
		c = fm6_comp(seq[i]);
		fm6_extend(e, &ik, ok.k, 0); // forward extension
		if (ok[c].x[2] != ik.x[2]) { // change of interval size
			int sc, ec;
			kv_push(fmintv_t, *curr, ik);
			if (i - x >= opt->min_l) {
				sc = ec2_score(opt, ok.k, c, &ec);
				if (sc - qual[i] > opt->for_qthres)
					c = ec, seq[i] = fm6_comp(c); // make a correction
			}
			if (ok[c].x[2] < opt->min_occ) break;
			// TODO: should I consider ambiguous bases more carefully?
		}
		ik = ok[c]; ik.info = i + 1;
	}

	if (i == len) kv_push(fmintv_t, *curr, ik);
	fm_reverse_fmivec(curr);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= 0; --i) { // backward search
		int sc = 0, ec;
		ec = c = seq[i];
		kv_resize(fmintv6_t, *tmp, prev->n);
		tmp->n = prev->n;
		for (j = 0; j < prev->n; ++j) // collect all the intervals at the current position
			fm6_extend(e, &prev->a[j], tmp->a[j].k, 1);
		for (j = 0; j < prev->n; ++j) { // check if we need to make a correction
			int sc2, ec2;
			sc2 = ec2_score(opt, tmp->a[j].k, c, &ec2);
		}
		if (sc - qual[i] > opt->rev_qthres)
			seq[i] = c = ec;
		for (j = 0; j < prev->n; ++j) { // update $curr
			uint64_t *ok = tmp->a[j].k;
			if (ok[c].x[2] >= opt->min_occ && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = p->info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	return ret;
}

