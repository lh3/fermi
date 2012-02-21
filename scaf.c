#include <zlib.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "priv.h"
#include "kstring.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_DECLARE(64, uint64_t, uint64_t)

#define A_THRES 20.
#define MIN_ISIZE 50

extern unsigned char seq_nt6_table[128];

typedef struct {
	int l, patched;
	double t;
	char *s;
} ext_t;

typedef struct {
	uint64_t k[2];
	ext_t ext[2];
	int len, nsr, maxo;
	uint8_t *seq;
	ku128_v reads;
	uint64_t dist[2], dist2[2];
	int64_t nei[2], nei2[2];
} utig_t;

typedef kvec_t(utig_t) utig_v;
typedef khash_t(64) hash64_t;

static utig_v *read_utig(const char *fn, int min_supp)
{
	int i, j, nsr;
	gzFile fp;
	kseq_t *kseq;
	utig_v *u;

	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	kseq = kseq_init(fp);
	u = calloc(1, sizeof(utig_v));
	while (kseq_read(kseq) >= 0) {
		char *q, *qq;
		long k[2];
		int beg, end;
		utig_t *p;

		if (kseq->comment.l == 0) continue; // no comments
		if ((q = strstr(kseq->comment.s, "UR:Z:")) == 0) continue; // no UR tag
		q += 5; // skip "UR:Z:"; jump to the first unmapped read (UR)
		qq = kseq->comment.s;
		nsr = strtol(qq, &qq, 10);
		if (nsr < min_supp) continue; // too few reads

		kv_pushp(utig_t, *u, &p);
		memset(p, 0, sizeof(utig_t));
		p->nei[0] = p->nei[1] = -1;
		sscanf(kseq->name.s, "%ld:%ld", &k[0], &k[1]);
		p->nsr = nsr;
		p->k[0] = k[0]; p->k[1] = k[1];
		beg = 0; end = kseq->seq.l;
		if (kseq->qual.l) { // trim unitigs covered by a single read
			for (i = 0; i < kseq->qual.l && kseq->qual.s[i] == 34; ++i);
			beg = i;
			for (i = kseq->qual.l - 1; i >= 0 && kseq->qual.s[i] == 34; --i);
			end = i + 1;
			if (beg >= end) beg = 0, end = kseq->seq.l;
		}
		p->len = end - beg;
		p->seq = calloc(1, end - beg + 1);
		strncpy((char*)p->seq, kseq->seq.s + beg, end - beg);
		for (i = 0; i < p->len; ++i)
			p->seq[i] = seq_nt6_table[(int)p->seq[i]];

		for (j = p->maxo = 0; j < 2; ++j) {
			if (*qq != '.') {
				while (isdigit(*qq) || *qq == '-') { // parse the neighbors
					int o;
					strtol(qq, &qq, 10); ++qq; // read rank
					o = strtol(qq, &qq, 10); ++qq;
					p->maxo = p->maxo > o? p->maxo : o;
				}
				++qq; // skip the tailing blank
			} else qq += 2;
		}

		while (isdigit(*q)) { // read mapping
			ku128_t x;
			x.x = strtol(q, &q, 10); ++q;
			x.y = (uint64_t)strtol(q, &q, 10)<<32; ++q;
			x.y |= strtol(q, &q, 10);
			kv_push(ku128_t, p->reads, x);
			if (*q++ == 0) break;
		}
	}
	kseq_destroy(kseq);
	gzclose(fp);
	return u;
}

static double cal_rdist(const utig_v *v)
{
	int j;
	uint64_t *srt;
	double rdist = -1.;
	int64_t i, sum_n_all, sum_n, sum_l;

	srt = calloc(v->n, 8);
	for (i = 0, sum_n_all = 0; i < v->n; ++i) {
		srt[i] = (uint64_t)v->a[i].nsr<<32 | i;
		sum_n_all += v->a[i].nsr;
	}
	ks_introsort_uint64_t(v->n, srt);
	for (j = 0; j < 2; ++j) {
		sum_n = sum_l = 0;
		for (i = v->n - 1; i >= 0; --i) {
			const utig_t *p = &v->a[srt[i]<<32>>32];
			if (rdist > 0. && (p->len - p->maxo) / rdist - p->nsr * M_LN2 < A_THRES) continue;
			sum_n += p->nsr;
			sum_l += p->len - p->maxo;
			if (sum_n >= sum_n_all * 0.5) break;
		}
		rdist = (double)sum_l / sum_n;
	}
	free(srt);
	return rdist;
}

static hash64_t *collect_nei(utig_v *v, int max_dist)
{
	int i, j, a, is_absent;
	hash64_t *h, *t;
	khint_t k;

	h = kh_init(64);
	for (i = 0; i < v->n; ++i) {
		utig_t *p = &v->a[i];
		int dist;
		for (j = 0; j < p->reads.n; ++j) {
			uint64_t idd = i<<1 | ((p->reads.a[j].x&1)^1);
			if (p->reads.a[j].x&1) dist = p->reads.a[j].y<<32>>32;
			else dist = p->len - (p->reads.a[j].y>>32);
			if (dist > max_dist) continue; // skip this read
			k = kh_put(64, h, p->reads.a[j].x>>1, &is_absent);
			if (is_absent) kh_val(h, k) = idd<<32 | dist;
			else kh_val(h, k) = 0; // mark delete
		}
	}
	for (k = 0; k != kh_end(h); ++k)
		if (kh_exist(h, k) && kh_val(h, k) == 0) kh_del(64, h, k); // now delete those that are marked "delete"

	t = kh_init(64);
	for (i = 0; i < v->n; ++i) {
		utig_t *q, *p = &v->a[i];
		for (a = 0; a < 2; ++a) {
			if (kh_n_buckets(t) >= 32) { // if t is too large, reallocate it
				kh_destroy(64, t);
				t = kh_init(64);
			} else kh_clear(64, t);
			for (j = 0; j < p->reads.n; ++j) {
				int dist;
				k = kh_get(64, h, p->reads.a[j].x>>1); // lookup the read
				if (k == kh_end(h) || (kh_val(h, k)>>32&1) != a) continue; // deleted or not in the right direction
				dist = (int32_t)kh_val(h, k);
				k = kh_get(64, h, p->reads.a[j].x>>1^1); // lookup the mate
				if (k == kh_end(h)) continue; // mate absent or deleted
				q = &v->a[kh_val(h, k)>>33];
				if (p == q) continue; // don't know how to deal with this case
				dist += (int32_t)kh_val(h, k);
				k = kh_put(64, t, kh_val(h, k)>>32, &is_absent);
				if (is_absent) kh_val(t, k) = 1ULL<<40 | dist;
				else kh_val(t, k) += 1ULL<<40 | dist;
			}
			for (k = 0; k != kh_end(t); ++k) { // write p->dist[a] and p->nei[a]
				if (!kh_exist(t, k) || kh_val(t, k)>>40 < 2) continue;
				if (kh_val(t, k) >= p->dist[a])
					p->dist2[a] = p->dist[a], p->nei2[a] = p->nei[a], p->dist[a] = kh_val(t, k), p->nei[a] = kh_key(t, k);
				else if (kh_val(t, k) >= p->dist2[a]) p->dist2[a] = kh_val(t, k), p->nei2[a] = kh_key(t, k);
			}
		}
	}
	kh_destroy(64, t);

	for (i = 0; i < v->n; ++i) { // test reciprocal best
		utig_t *q, *p = &v->a[i];
		for (a = 0; a < 2; ++a) {
			if (p->nei[a] < 0) continue;
			q = &v->a[p->nei[a]>>1];
			if (q->nei[p->nei[a]&1] != (i<<1|a))
				p->dist[a] = 0; // we should not set p->nei[a]=-1 at this time, as it may interfere with others
		}
	}
	for (i = 0; i < v->n; ++i) {
		utig_t *p = &v->a[i];
		if (p->dist[0] == 0) p->nei[0] = -2;
		if (p->dist[1] == 0) p->nei[1] = -2;
	}
	for (i = 0; i < v->n; ++i) { // update nei2[] as nei[] is now changed
		utig_t *p = &v->a[i];
		for (a = 0; a < 2; ++a) {
			if (p->nei[a] >= 0 && p->nei2[a] >= 0) {
				utig_t *q = &v->a[p->nei2[a]>>1];
				if (q->nei[p->nei2[a]&1] < 0) p->nei2[a] = -3;
			}
		}
	}
	return h;
}

#if 1
static void debug_utig(utig_v *v, uint32_t idd, double rdist)
{
	int a = idd&1;
	utig_t *p = &v->a[idd>>1];
	if (p->nei[a] >= 0) {
		fprintf(stderr, "%d[%ld:%ld]\t%d:%d:%f\t%ld\t%d:%ld\t%ld\t%d:%ld\t%d\t%d\t%g\n", idd, (long)p->k[0], (long)p->k[1], p->len, p->nsr, (p->len - p->maxo)/rdist - M_LN2 * p->nsr,
				(long)p->nei[a], (int)(p->dist[a]>>40), (long)((double)(p->dist[a]<<24>>24)/(p->dist[a]>>40) + .499),
				(long)p->nei2[a], (int)(p->dist2[a]>>40), (long)((double)(p->dist2[a]<<24>>24)/(p->dist2[a]>>40) + .499),
				p->ext[a].patched, p->ext[a].l, p->ext[a].t);
	}
}
#endif

/***************************************
 * Gamma and incomplete Beta functions *
 ***************************************/

#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

double kf_lgamma(double z)
{
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

static double kf_betai_aux(double a, double b, double x)
{
	double C, D, f;
	int j;
	if (x == 0.) return 0.;
	if (x == 1.) return 1.;
	f = 1.; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	for (j = 1; j < 200; ++j) {
		double aa, d;
		int m = j>>1;
		aa = (j&1)? -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))
			: m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m));
		D = 1. + aa * D;
		if (D < KF_TINY) D = KF_TINY;
		C = 1. + aa / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(kf_lgamma(a+b) - kf_lgamma(a) - kf_lgamma(b) + a * log(x) + b * log(1.-x)) / a / f;
}
double kf_betai(double a, double b, double x)
{
	return x < (a + 1.) / (a + b + 2.)? kf_betai_aux(a, b, x) : 1. - kf_betai_aux(b, a, 1. - x);
}

/***************
 * Gap closure *
 ***************/

static inline void end_seq(kstring_t *str, const utig_t *p, int is3, int is_2nd, int max_dist)
{
	int ori_l = str->l;
	if (p->len > max_dist) {
		if (is3) kputsn((char*)p->seq + (p->len - max_dist), max_dist, str);
		else kputsn((char*)p->seq, max_dist, str);
	} else kputsn((char*)p->seq, p->len, str);
	if ((!is3) ^ (!!is_2nd)) seq_revcomp6(str->l - ori_l, (uint8_t*)str->s + ori_l);
	kputc(0, str);
}

static int add_seq(const rld_t *e, const hash64_t *h, const utig_t *p, int64_t idd, kstring_t *str, kstring_t *tmp)
{
	int j, max_len;
	for (j = max_len = 0; j < p->reads.n; ++j) {
		khint_t k = kh_get(64, h, p->reads.a[j].x>>1^1);
		if (k != kh_end(h) && (idd < 0 || kh_val(h, k)>>32 == idd)) {
			assert(p->reads.a[j].x < e->mcnt[1]);
			fm_retrieve(e, p->reads.a[j].x, tmp);
			if (tmp->l > max_len) max_len = tmp->l;
			seq_reverse(tmp->l, (uint8_t*)tmp->s);
			kputsn(tmp->s, tmp->l + 1, str);
		}
	}
	return max_len;
}

static void print_seq(char *s)
{
	for (; *s; ++s) putchar("$ACGTN"[(int)*s]);
	putchar('\n');
}

static double compute_t(const hash64_t *h, const utig_v *v, uint32_t idd, int l, double mu)
{
	utig_t *p = &v->a[idd>>1];
	int j, n, dist;
	int64_t sum, sum2;
	double t, avg;
	if (p->nei[idd&1] < 0) return 0.0;
	sum = sum2 = 0; n = 0;
	for (j = 0; j < p->reads.n; ++j) {
		khint_t k = kh_get(64, h, p->reads.a[j].x>>1);
		if (k == kh_end(h)) continue;
		dist = kh_val(h, k)<<32>>32;
		k = kh_get(64, h, p->reads.a[j].x>>1^1);
		if (k == kh_end(h) || kh_val(h, k)>>32 != p->nei[idd&1]) continue;
		dist += kh_val(h, k)<<32>>32;
		dist += l;
		++n; sum += dist; sum2 += dist * dist;
	}
	assert(n >= 2);
	avg = (double)sum / n;
	t = sqrt(((double)sum2 / n - avg * avg) / (n - 1)); // std.dev. / sqrt(n)
	t = (avg - mu) / t; // student's t
	--n; // n is now the degree of freedom
	if (n > 50) n = 50; // avoid a too stringent P-value
	return kf_betai(.5*n, .5, n/(n+t*t));
}

static ext_t assemble(int l, char *s, int max_len, char *const t[2])
{
	mag_t *g;
	magv_t *p;
	int j, max_j, old_verbose = fm_verbose;
	char *q, *r;
	ext_t e;

	fm_verbose = 1;
	memset(&e, 0, sizeof(ext_t));
	g = fm6_api_unitig(max_len/3. < 17? max_len/3. : 17, l, s);
	mag_g_rm_vext(g, max_len + 1, 100);
	mag_g_merge(g, 0);
	mag_g_simplify_bubble(g, 25, max_len * 2);
	mag_g_pop_simple(g, 10., 0.15, 1); // FIXME: always agressive?
	for (j = max_len = 0; j < g->v.n; ++j)
		if (g->v.a[j].len > max_len)
			max_len = g->v.a[j].len, max_j = j;
	p = &g->v.a[max_j];
//	print_seq(t[0]); print_seq(t[1]); mag_g_print(g);
	q = strstr(p->seq, t[0]);
	if (q == 0) {
		seq_revcomp6(p->len, (uint8_t*)p->seq);
		q = strstr(p->seq, t[0]);
	}
	if (q) {
		if ((r = strstr(p->seq, t[1])) > q) { // gap patched
			int tmp = strlen(t[0]);
			e.patched = 1;
			e.l = r - (q + tmp);
			if (e.l > 0) {
				e.s = calloc(1, e.l + 1);
				strncpy(e.s, p->seq + tmp, e.l);
			}
		}
	}
	mag_g_destroy(g);
	fm_verbose = old_verbose;
	return e;
}

static void patch_gap(const rld_t *e, const hash64_t *h, utig_v *v, uint32_t iddp, int max_dist, double avg)
{
	uint32_t iddq;
	utig_t *p, *q;
	kstring_t str, rd;
	int max_len, pl, i;
	char *t[2];
	ext_t ext;

	p = &v->a[iddp>>1];
	if (p->nei[iddp&1] < 0) return; // no neighbor
	iddq = p->nei[iddp&1];
	if (iddp > iddq) return; // avoid doing local assembly twice
	q = &v->a[iddq>>1];
	if (q->nei[iddq&1] != iddp) return; // not reciprocal best

	str.s = rd.s = 0; str.m = rd.m = 0;
	for (i = 0; i < 2; ++i) {
		str.l = rd.l = 0;
		end_seq(&str, p, iddp&1, 0, max_dist); pl = str.l;
		end_seq(&str, q, iddq&1, 1, max_dist);
		max_len = add_seq(e, h, p, i? -1 : iddq, &str, &rd); // the first round, using reads from unitigs only
		add_seq(e, h, q, i? -1 : iddp, &str, &rd); // the second round, using all unpaired reads
		t[0] = str.s; t[1] = str.s + pl;
		ext = assemble(str.l, str.s, max_len, t);
		if (ext.patched) {
			ext.t = compute_t(h, v, iddp, ext.l, avg);
			if (ext.t > 1e-6) {
				p->ext[iddp&1] = q->ext[iddq&1] = ext;
				break;
			}
		}
	}

	free(str.s); free(rd.s);
}

/**********
 * Portal *
 **********/

void mag_scaf_core(const rld_t *e, const char *fn, double avg, double std, int min_supp)
{
	utig_v *v;
	double rdist;
	hash64_t *h;
	int i, max_dist;

	max_dist = (int)(avg + 2. * std + .499);
	v = read_utig(fn, min_supp);
	rdist = cal_rdist(v);
	h = collect_nei(v, max_dist);
//	patch_gap(e, h, v, 64, max_dist); debug_utig(v, 64, rdist); return;
	for (i = 0; i < v->n; ++i) {
		patch_gap(e, h, v, i<<1|0, max_dist, avg);
		patch_gap(e, h, v, i<<1|1, max_dist, avg);
	}
	for (i = 0; i < v->n; ++i) {
		debug_utig(v, i<<1|0, rdist);
		debug_utig(v, i<<1|1, rdist);
	}
	kh_destroy(64, h);
}
