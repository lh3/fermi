#include "priv.h"
#include "kvec.h"

void fmc_opt_init(fmec2opt_t *opt)
{
	opt->min_l = 19;
	opt->min_occ_f = 6;
	opt->min_occ_r = 3;

	opt->max_pen = 60;
	opt->max_d = 7;
	opt->diff_factor = 13;
	opt->ratio_factor = 10.;
}

static inline void ec_cns_gen(const fmec2opt_t *opt, const fmintv_t k[6], uint8_t q[4])
{
	int c;
	int64_t sum = 0;
	for (c = 1; c <= 4; ++c) sum += k[c].x[2];
	for (c = 1; c <= 4; ++c) {
		int64_t d0, d1;
		d0 = k[c].x[2]; d1 = sum - d0;
		if (d0 < d1) {
			int p, r;
			p = ((d1 < opt->max_d? d1 : opt->max_d) - d0) * opt->diff_factor;
			p = p > 1? p : 1;
			p = p < opt->max_pen? p : opt->max_pen;
			r = d0? (int)(opt->ratio_factor * d1 / d0 + .499) : opt->max_pen;
			q[c-1] = p < r? p : r;
		} else q[c-1] = 0;
	}
}

static inline void ec_cns_upd(uint8_t q[4], const uint8_t q1[4], int is_comp)
{
	int c;
	if (is_comp) {
		for (c = 0; c < 4; ++c)
			q[c] = q1[3-c];
	} else {
		for (c = 0; c < 4; ++c)
			q[c] = q1[c];
	}
}

static inline int ec_cns_call(const uint8_t q[4], int base, int qual)
{
	int min, min_c, min2, c, iq[4];
	for (c = 0; c < 4; ++c) iq[c] = q[c];
	if (base >= 0 && base <= 3) iq[base] -= qual;
	for (c = 0, min_c = -1, min = min2 = 255; c < 4; ++c)
		if (min > iq[c]) min2 = min, min = iq[c], min_c = c;
		else if (min2 > iq[c]) min2 = iq[c];
	return min_c << 8 | (min2 - min);
}

typedef struct {
	uint8_t pen[4];
	uint8_t ob, oq, eb, eq;
} ecseq_t;

static ecseq_t *ec_seq_gen(int l_seq, const char *seq, const char *qual)
{
	int i;
	ecseq_t *s;
	s = calloc(l_seq, sizeof(ecseq_t));
	for (i = 0; i < l_seq; ++i) {
		uint32_t *p = (uint32_t*)s[i].pen;
		*p = 0xffffffffU;
		s[i].ob = seq_nt6_table[(int)seq[i]];
		s[i].eq = s[i].oq = qual[i] - 33;
	}
	return s;
}

static void ec_seq_rev(int l_seq, ecseq_t *seq)
{
	int i;
	for (i = 0; i < l_seq>>1; ++i) {
		ecseq_t tmp, *s = &seq[i], *t = &seq[l_seq - 1 - i];
		tmp = *s; *s = *t; *t = tmp;
	}
	for (i = 0; i < l_seq; ++i) {
		ecseq_t *si = &seq[i];
		uint32_t *p = (uint32_t*)si->pen;
		si->ob = fm6_comp(si->ob); si->eb = fm6_comp(si->eb);
		*p = (*p & 0x0000FFFFU) << 16 | *p >> 16;
		*p = (*p & 0x00FF00FFU) << 8  | (*p & 0xFF00FF00U) >> 8;
	}
}

static int ec2_best_intv(const fmec2opt_t *opt, int i, int n, fmintv_t *a)
{
	int k, max_k;
	int64_t sum = 0, s;
	double max, x;
	for (k = n - 1; k >= 0; --k)
		if (a[k].info - i >= opt->min_l) break;
	n = k + 1;
	if (n == 0) return -1;
	if (n == 1) return 0;
	// 2-mean
	for (k = 0; k < n; ++k) sum += a[k].x[2];
	for (k = 1, s = 0, max = 0., max_k = 0; k < n; ++k) {
		s += a[k-1].x[2];
		x = (double)s * s / k + (double)(sum - s) * (sum - s) / (n - k);
		if (max < x) max = x, max_k = k;
	}
	return max_k - 1;
}

static int ec2_core(const fmec2opt_t *opt, const rld_t *e, int l_seq, ecseq_t *seq, int x, fmintv_v *prev, fmintv_v *curr, fmintv6_v *tmp)
{ // this function is similar to fm6_smem1_core()
	int i, j, ret;
	fmintv_t ik;
	fmintv_v *swap;
	fmintv6_t ok;

	if (fm_verbose >= 5) fprintf(stderr, "x=%d\n", x);
	fm6_set_intv(e, seq[x].ob, ik);
	ik.info = x + 1;
	for (i = x + 1, curr->n = 0; i < l_seq; ++i) { // forward search
		ecseq_t *si = &seq[i];
		int c;
		fm6_extend(e, &ik, ok.k, 0); // forward extension
		if (0&&i - x >= opt->min_l) { // then check the multi-alignment
			int call;
			uint8_t q[4];
			ec_cns_gen(opt, ok.k, q);
			ec_cns_upd(si->pen, q, 1);
			call = ec_cns_call(si->pen, fm6_comp(si->ob) - 1, si->oq);
			si->eb = (call >> 8) + 1;
			si->eq = call & 0xff;
			if (fm_verbose >= 5 && si->eb != si->ob)
				fprintf(stderr, "[F,%d,%c]\tpen[4]={%d,%d,%d,%d}\teb=%c\teq=%d\n", i, "$ACGTN"[si->ob], q[0], q[1], q[2], q[3], "ACGT"[call>>8], call&0xff);
		}
		c = si->eb? si->eb : si->ob;
		c = fm6_comp(c);
		if (ok.k[c].x[2] != ik.x[2]) { // change of interval size
			kv_push(fmintv_t, *curr, ik);
			if (ok.k[c].x[2] < opt->min_occ_f) break;
		}
		ik = ok.k[c]; ik.info = i + 1;
	}

	if (i == l_seq) kv_push(fmintv_t, *curr, ik);
	fm_reverse_fmivec(curr);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= 0; --i) { // backward search
		ecseq_t *si = &seq[i];
		int c, best_intv;
		kv_resize(fmintv6_t, *tmp, prev->n);
		tmp->n = prev->n;
		for (j = 0; j < prev->n; ++j) // collect all the intervals at the current position
			fm6_extend(e, &prev->a[j], tmp->a[j].k, 1);
		c = si->eb? si->eb : si->ob;
		best_intv = ec2_best_intv(opt, i, curr->n, curr->a);
		if (best_intv >= 0) {
			uint8_t q[4];
			int call;
			ec_cns_gen(opt, tmp->a[best_intv].k, q);
			ec_cns_upd(si->pen, q, 0);
			call = ec_cns_call(si->pen, si->ob - 1, si->oq);
			si->eb = c = (call >> 8) + 1;
			si->eq = call & 0xff;
			//if (fm_verbose >= 5 && si->eb != si->ob)
			if (fm_verbose >= 5)
				fprintf(stderr, "[R,%d,%c]\tpen[4]={%d,%d,%d,%d}\teb=%c\teq=%d\tbest=%d\n", i, "$ACGTN"[si->ob], si->pen[0], si->pen[1], si->pen[2], si->pen[3], "ACGT"[call>>8], call&0xff, best_intv);
		}
		for (j = 0, curr->n = 0; j < prev->n; ++j) { // update $curr
			fmintv_t *ok = tmp->a[j].k;
			if (ok[c].x[2] >= opt->min_occ_r && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = prev->a[j].info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	//fprintf(stderr, "ret=%d\n", ret);

	return ret;
}

int fm6_ec2_core(const fmec2opt_t *opt, const rld_t *e, int l_seq, char *seq, char *qual, fmintv_v *prev, fmintv_v *curr, fmintv6_v *tmp)
{
	int x;
	ecseq_t *s;
	s = ec_seq_gen(l_seq, seq, qual);
	// error correction
	x = 0;
	while ((x = ec2_core(opt, e, l_seq, s, x, prev, curr, tmp)) < l_seq);
	ec_seq_rev(l_seq, s);
	x = 0;
	while ((x = ec2_core(opt, e, l_seq, s, x, prev, curr, tmp)) < l_seq);
	ec_seq_rev(l_seq, s);
	// modify $seq and $qual
	for (x = 0; x < l_seq; ++x) {
		ecseq_t *sx = &s[x];
		seq[x] = sx->eb == 0 || sx->ob == sx->eb? "$ACGTN"[sx->ob] : "$acgtn"[sx->eb];
		qual[x] = sx->oq + 33;
	}
	free(s);
	return 0;
}
