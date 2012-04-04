#include <zlib.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "priv.h"
#include "kstring.h"
#include "kvec.h"
#include "ksw.h"
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
	double A;
	int len, nsr, maxo;
	uint16_t deleted, excluded;
	uint8_t *seq;
	ku128_v reads;
	uint64_t dist[2], dist2[2];
	int64_t nei[2], nei2[2];
} utig_t;

typedef kvec_t(utig_t) utig_v;
typedef khash_t(64) hash64_t;

/*************
 * Basic I/O *
 *************/

static utig_v *read_utig(const char *fn)
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

		kv_pushp(utig_t, *u, &p);
		memset(p, 0, sizeof(utig_t));
		p->nei[0] = p->nei[1] = p->nei2[0] = p->nei2[1] = -1;
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
			int b, e;
			x.x = strtol(q, &q, 10); ++q;
			b = strtol(q, &q, 10); ++q;
			e = strtol(q, &q, 10);
			x.y = (uint64_t)(b > beg? b - beg : 0)<<32 | (e - beg < p->len? e - beg : p->len);
			kv_push(ku128_t, p->reads, x);
			if (*q++ == 0) break;
		}
	}
	kseq_destroy(kseq);
	gzclose(fp);
	return u;
}

static void utig_destroy(utig_v *v)
{
	int i, a;
	for (i = 0; i < v->n; ++i) {
		utig_t *p = &v->a[i];
		free(p->seq); free(p->reads.a);
		for (a = 0; a < 2; ++a)
			if (i<<1 < p->nei[a]) free(p->ext[a].s);
	}
	free(v->a); free(v);
}

static void debug_utig(utig_v *v, uint32_t idd)
{
	int b, a = idd&1;
	utig_t *q, *p = &v->a[idd>>1];
	fprintf(stderr, "LK\t%u:%d\t%ld\t%d\t%d\t%.2f", idd>>1, idd&1, (long)p->k[a], p->len, p->nsr, p->A);
	if (p->nei[a] >= 0) {
		q = &v->a[p->nei[a]>>1];
		b = p->nei[a]&1;
		fprintf(stderr, "\t%ld\t%d:%d", (long)q->k[b], (int)(p->dist[a]>>40), (int)(p->dist[a]<<24>>24));
		fprintf(stderr, "\t%d:%d:%.1e", p->ext[a].patched, p->ext[a].l, p->ext[a].t);
	}
	if (p->nei2[a] >= 0) {
		q = &v->a[p->nei2[a]>>1];
		b = p->nei2[a]&1;
		fprintf(stderr, "\t%ld\t%d:%d", (long)q->k[b], (int)(p->dist2[a]>>40), (int)(p->dist2[a]<<24>>24));
	}
	fputc('\n', stderr);
}

/*********************************
 * Compute A and connect unitigs *
 *********************************/

static double cal_rdist(utig_v *v)
{
	int j, n_ovlp, avg_ovlp;
	uint64_t *srt;
	double rdist = -1.;
	int64_t i, sum_n_all, sum_n, sum_l, sum_ovlp;

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

	sum_ovlp = 0; n_ovlp = 0;
	for (i = 0; i < v->n; ++i)
		if (v->a[i].maxo) ++n_ovlp, sum_ovlp += v->a[i].maxo;
	avg_ovlp = (int)((double)sum_ovlp / n_ovlp + .499);
	for (i = 0; i < v->n; ++i) {
		utig_t *p = &v->a[i];
		p->A = (p->len - (p->maxo? p->maxo : avg_ovlp)) / rdist - p->nsr * M_LN2;
	}
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
		if (p->excluded) continue;
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
				if (!kh_exist(t, k) || kh_val(t, k)>>40 < 1) continue;
				if (kh_val(t, k) >= p->dist[a])
					p->dist2[a] = p->dist[a], p->nei2[a] = p->nei[a], p->dist[a] = kh_val(t, k), p->nei[a] = kh_key(t, k);
				else if (kh_val(t, k) >= p->dist2[a]) p->dist2[a] = kh_val(t, k), p->nei2[a] = kh_key(t, k);
			}
		}
	}
	kh_destroy(64, t);

	// change the lower 40-bit as the average distance, instead of sum;
	for (i = 0; i < v->n; ++i) {
		utig_t *p = &v->a[i];
		for (a = 0; a < 2; ++a) {
			if (p->dist[a]) p->dist[a] = p->dist[a]>>40<<40 | (int)((double)(p->dist[a]<<24>>24) / (p->dist[a]>>40) + .499);
			if (p->dist2[a]) p->dist2[a] = p->dist2[a]>>40<<40 | (int)((double)(p->dist2[a]<<24>>24) / (p->dist2[a]>>40) + .499);
		}
	}
	return h;
}

static void resolve_contained(utig_v *v, uint32_t id, double avg, double std, int pr_link) // FIXME: only works in simple cases; topological sorting is better
{
	utig_t *p = &v->a[id], *q[2];
	int d_long, d_short, a;
	if (p->excluded || p->nei[0] < 0 || p->nei[1] < 0 || p->nei2[0] >= 0 || p->nei2[1] >= 0) return;
	q[0] = &v->a[p->nei[0]>>1]; q[1] = &v->a[p->nei[1]>>1];
	if (q[0]->nei2[p->nei[0]&1] < 0 || q[1]->nei2[p->nei[1]&1] < 0) return;
	if (q[1]->nei[p->nei[1]&1] != p->nei[0] && q[1]->nei2[p->nei[1]&1] != p->nei[0]) return;
	if (q[0]->nei[p->nei[0]&1] == p->nei[1]) {
		d_long = (int)(avg - (q[0]->dist[p->nei[0]&1]<<24>>24) + .499);
	} else if (q[0]->nei2[p->nei[0]&1] == p->nei[1]) {
		d_long = (int)(avg - (q[0]->dist2[p->nei[0]&1]<<24>>24) + .499);
	} else return;
	d_short = (int)(2*avg - (p->dist[0]<<24>>24) - (p->dist[1]<<24>>24) + p->len + .499);
	if (abs(d_long - d_short) < std) {
		if (pr_link) {
			fprintf(stderr, "CT\t%ld:%ld\t%d\t%d\n", (long)p->k[0], (long)p->k[1], d_long, d_short);
			// break the link between q[0] and q[1]
			for (a = 0; a < 2; ++a) {
				if (q[a]->nei[p->nei[a]&1] == p->nei[a^1]) {
					q[a]->nei[p->nei[a]&1] = q[a]->nei2[p->nei[a]&1];
					q[a]->dist[p->nei[a]&1] = q[a]->dist2[p->nei[a]&1];
				}
				q[a]->nei2[p->nei[a]&1] = -4;
				q[a]->dist2[p->nei[a]&1] = 0;
			}
		}
	}
}

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

static int add_seq(const rld_t *e, const hash64_t *h, const utig_t *p, kstring_t *str, kstring_t *tmp, int64_t idd_self, int64_t idd_mate)
{
	int j, max_len;
	for (j = max_len = 0; j < p->reads.n; ++j) {
		khint_t k = kh_get(64, h, p->reads.a[j].x>>1);
		if (k == kh_end(h) || kh_val(h, k)>>32 != idd_self) continue; // make sure: 1) not duplicated; 2) in the right direction; 3) close to the end
		if (idd_mate >= 0) {
			k = kh_get(64, h, p->reads.a[j].x>>1^1);
			if (k == kh_end(h) || kh_val(h, k)>>32 != idd_mate) continue;
		}
		assert((p->reads.a[j].x^3) < e->mcnt[1]);
		fm_retrieve(e, p->reads.a[j].x^3, tmp); // retrieve the mate
		if (tmp->l > max_len) max_len = tmp->l;
		seq_reverse(tmp->l, (uint8_t*)tmp->s); // sequence returned by fm_retrieve() are reversed (but not complemeted)
		kputsn(tmp->s, tmp->l + 1, str);
	}
	return max_len;
}

static double correct_mean(double l, double mu, double sigma)
{
	double x, y, z;
	x = (l - mu) / sigma;
	y = M_SQRT2 / M_2_SQRTPI * erfc(x * M_SQRT1_2);
	z = exp(-.5 * x * x);
	return mu + sigma * y / (z - x * y);
}

static double compute_t(const hash64_t *h, const utig_v *v, uint32_t idd, int l, double mu, double sigma, int max_len)
{
	utig_t *p = &v->a[idd>>1];
	int j, n, dist;
	int64_t sum, sum2;
	double t, avg, mu_;
	if (p->nei[idd&1] < 0) return 0.0;
	sum = sum2 = 0; n = 0;
	mu_ = correct_mean(2 * max_len + l, mu, sigma);
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
	t = (avg - mu_) / t; // student's t
	--n; // n is now the degree of freedom
	if (n > 50) n = 50; // avoid a too stringent P-value
	return kf_betai(.5*n, .5, n/(n+t*t));
}

static ext_t assemble(int l, char *s, int max_len, char *const t[2])
{
	mag_t *g;
	magv_t *p;
	int j, max_j;
	char *q, *r;
	ext_t e;

	memset(&e, 0, sizeof(ext_t));
	//printf(">0\n");for(j=0;j<l-1;++j)if(s[j]==0)printf("\n>%d\n",j);else putchar("$ACGTN"[(int)s[j]]);putchar('\n');exit(0);
	g = fm6_api_unitig(max_len/3. < 17? max_len/3. : 17, l, s);
	mag_g_merge(g, 1); // FIXME: this to remove multi-edges, which is likely to introduce small scale errors...
	mag_g_rm_vext(g, max_len * 1.1, 4);
	mag_g_simplify_bubble(g, 25, max_len * 2);
	mag_g_pop_simple(g, 10., 0.15, 1); // FIXME: always agressive?
	mag_g_rm_edge(g, 0, 0.8, max_len * 1.1, 5);
	mag_g_merge(g, 1);
	mag_g_rm_vext(g, max_len * 1.1, 100);
	mag_g_merge(g, 0);
	mag_g_simplify_bubble(g, 25, max_len * 2);
	mag_g_pop_simple(g, 10., 0.15, 1); // FIXME: always agressive?
	for (j = max_len = 0, max_j = -1; j < g->v.n; ++j)
		if (g->v.a[j].len > max_len)
			max_len = g->v.a[j].len, max_j = j;
	if (max_j >= 0) { // sometimes the whole graph can be empty
		p = &g->v.a[max_j];
		// mag_g_print(g);
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
	}
	mag_g_destroy(g);
	return e;
}

#define MAX_DROP 7
#define SCORE_THRES 13

static void patch_gap(const rld_t *e, const hash64_t *h, utig_v *v, uint32_t iddp, int min_supp, int max_dist, double avg, double std)
{
	uint32_t iddq;
	utig_t *p, *q;
	kstring_t str, rd;
	int max_len, pl, ql, i, dist1, dist2;
	char *t[2];
	ext_t ext;

	p = &v->a[iddp>>1];
	if (p->nei[iddp&1] < 0 || p->dist[iddp&1]>>40 < min_supp) return; // no neighbor
	iddq = p->nei[iddp&1];
	if (iddp >= iddq) return; // avoid doing local assembly twice
	q = &v->a[iddq>>1];
	if (q->nei[iddq&1] != iddp) return; // not reciprocal best

	dist1 = p->dist[iddp&1]>>40; dist2 = 0;
	if (p->nei2[iddp&1] >= 0) dist2 = p->dist2[iddp&1]>>40;
	if (q->nei2[iddq&1] >= 0) dist2 = dist2 > q->dist2[iddq&1]>>40? dist2 : q->dist2[iddq&1]>>40;
	if (dist2 >= min_supp || (double)dist2 / dist1 >= 1./min_supp) return;

	str.s = rd.s = 0; str.m = rd.m = 0;
	for (i = 0; i < 2; ++i) {
		str.l = rd.l = 0;
		end_seq(&str, p, iddp&1, 0, max_dist); pl = str.l;
		end_seq(&str, q, iddq&1, 1, max_dist); ql = str.l - pl;
		max_len = add_seq(e, h, p, &str, &rd, iddp, i? -1L : (int64_t)iddq);
		add_seq(e, h, q, &str, &rd, iddq, i? -1L : (int64_t)iddp);
		t[0] = str.s; t[1] = str.s + pl;
		ext = assemble(str.l, str.s, max_len, t);
		if (ext.patched && ext.l + p->len > 0 && ext.l + q->len > 0) {
			ext.t = compute_t(h, v, iddp, ext.l, avg, std, max_len);
			if (i == 0 && ext.t > 1e-5) {
				p->ext[iddp&1] = q->ext[iddq&1] = ext;
				break;
			} else if (i == 1 && ext.t > 1e-10) p->ext[iddp&1] = q->ext[iddq&1] = ext;
		}
	}
	if (ext.patched == 0 && p->dist[iddp&1]<<24>>24 > avg) { // another try in case there are heterozygotes in the overlap
		int j, k, drop[2], max_drop, min_drop;
		int8_t mat[25];
		kswr_t a;
		for (i = k = 0; i < 5; ++i)
			for (j = 0; j < 5; ++j)
				mat[k++] = i == j? 1 : -3;
		a = ksw_align(ql - 1, (uint8_t*)t[1], pl - 1, (uint8_t*)t[0], 5, mat, 5, 2, KSW_XSTART, 0);
		drop[0] = a.qb; drop[1] = (pl - 1) - (a.te + 1);
		max_drop = drop[0] > drop[1]? drop[0] : drop[1];
		min_drop = drop[0] < drop[1]? drop[0] : drop[1];
		if (min_drop == 0 && max_drop < MAX_DROP && a.score >= SCORE_THRES + max_drop) { // an end-to-end alignment
			int lp = a.te + 1 - a.tb + drop[0] + drop[1];
			int lq = a.qe + 1 + drop[0] + drop[1];
			if (lp < p->len && lq < q->len) {
				p->ext[iddp&1].l = -lp;
				q->ext[iddq&1].l = -lq;
				p->ext[iddp&1].patched = q->ext[iddq&1].patched = 1;
				p->ext[iddp&1].t = q->ext[iddq&1].t = compute_t(h, v, iddp, p->ext[iddp&1].l, avg, std, max_len);
			}
		}
		if (!p->ext[iddp&1].patched) fprintf(stderr, "SW\t%ld\t%ld\t%d\t%d\t%d\n", (long)p->k[iddp&1], (long)q->k[iddq&1], drop[0], drop[1], a.score);
		//fprintf(stderr, "%c, %d, (%d, %d, %d), (%d, %d, %d)\n", "NY"[p->ext[iddp&1].patched], a.score, ql, a.qb, a.qe+1, pl-1, a.tb, a.te+1);
	}
	free(str.s); free(rd.s);
}

/****************
 * Join unitigs *
 ****************/

static void find_path1(utig_v *v, ku64_v *path, double a_thres, double p_thres)
{
	if (path->n == 0) return;
	for (;;) {
		uint32_t iddq, idd = path->a[path->n - 1]; // the last idd
		utig_t *q, *p = &v->a[idd>>1];
		if (p->nei[idd&1] < 0 || p->ext[idd&1].patched == 0 || p->ext[idd&1].t < p_thres) break;
		iddq = p->nei[idd&1];
		q = &v->a[iddq>>1];
		if (q->deleted || q->A < a_thres) break;
		kv_push(uint64_t, *path, iddq);
		kv_push(uint64_t, *path, iddq^1);
		q->deleted = 1;
	}
}

static void find_path(utig_v *v, uint32_t id, ku64_v *path, double a_thres, double p_thres)
{
	utig_t *p = &v->a[id];
	path->n = 0;
	if (p->deleted) return; // already used in other paths
	kv_push(uint64_t, *path, id<<1|0);
	kv_push(uint64_t, *path, id<<1|1);
	p->deleted = 1;
	if (p->A >= a_thres) {
		int i;
		find_path1(v, path, a_thres, p_thres);
		for (i = 0; i < path->n>>1; ++i) {
			uint64_t tmp;
			tmp = path->a[i];
			path->a[i] = path->a[path->n - 1 - i];
			path->a[path->n - 1 - i] = tmp;
		}
		find_path1(v, path, a_thres, p_thres);
	}
}

static void make_scaftigs(utig_v *v, double a_thres, double p_thres)
{
	int i, j;
	ku64_v path;
	kstring_t ctg;
	kv_init(path);
	ctg.l = ctg.m = 0; ctg.s = 0;
	for (i = 0; i < v->n; ++i) {
		find_path(v, i, &path, a_thres, p_thres);
		if (path.n) {
			int nsr = 0;
			utig_t *beg, *end;
			ctg.l = 0;
			assert(path.n % 2 == 0);
			for (j = 0; j < path.n; j += 2) {
				uint32_t idd = path.a[j], ndir = (idd&1)^1;
				int ori_l = ctg.l;
				utig_t *p = &v->a[idd>>1];
				nsr += p->nsr;
				kputsn((char*)p->seq, p->len, &ctg);
				if (idd&1) seq_revcomp6(ctg.l - ori_l, (uint8_t*)ctg.s + ori_l);
				if (j == path.n - 2) break;
				assert(p->ext[ndir].patched);
				if (p->ext[ndir].l > 0) {
					ori_l = ctg.l;
					kputsn(p->ext[ndir].s, p->ext[ndir].l, &ctg);
					if (path.a[j+2] < path.a[j])
						seq_revcomp6(ctg.l - ori_l, (uint8_t*)ctg.s + ori_l);
				} else ctg.l += p->ext[ndir].l;
			}
			for (j = 0; j < ctg.l; ++j)
				ctg.s[j] = "$ACGTN"[(int)ctg.s[j]];
			beg = &v->a[path.a[0]>>1]; end = &v->a[path.a[path.n-1]>>1];
			printf(">%ld:%ld\t%ld\t%d\t%.2f\n", (long)beg->k[path.a[0]&1], (long)end->k[path.a[path.n-1]&1], path.n/2, nsr, path.n > 2? 100.0 : beg->A);
			puts(ctg.s);
		}
	}
	free(path.a); free(ctg.s);
}

/*********************************
 * Multithreading local assembly *
 *********************************/

typedef struct {
	int start, step, max_dist;
	const rld_t *e;
	const hash64_t *h;
	const fmscafopt_t *opt;
	utig_v *v;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int64_t i;
	for (i = w->start; i < w->v->n; i += w->step) {
		patch_gap(w->e, w->h, w->v, i<<1|0, w->opt->min_supp, w->max_dist, w->opt->avg, w->opt->std);
		patch_gap(w->e, w->h, w->v, i<<1|1, w->opt->min_supp, w->max_dist, w->opt->avg, w->opt->std);
	}
	return 0;
}

/**********
 * Portal *
 **********/

void mag_scaf_core(const rld_t *e, const char *fn, const fmscafopt_t *opt, int n_threads)
{
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	utig_v *v;
	double rdist, t, treal;
	hash64_t *h;
	int i, max_dist, old_verbose;

	max_dist = (int)(opt->avg + 2. * opt->std + .499);
	t = cputime();
	v = read_utig(fn);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] read unitigs in %.3f sec\n", __func__, cputime() - t);
	t = cputime();
	rdist = cal_rdist(v);
	for (i = 0; i < v->n; ++i)
		if (v->a[i].A < opt->a_thres) v->a[i].excluded = 1;
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] rdist = %.3f, computed in %.3f sec\n", __func__, rdist, cputime() - t);
	t = cputime();
	h = collect_nei(v, max_dist);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] paired unitigs in %.3f sec\n", __func__, cputime() - t);
	for (i = 0; i < v->n; ++i)
		resolve_contained(v, i, opt->avg, opt->std, opt->pr_links);

//	patch_gap(e, h, v, 296, opt->min_supp, max_dist, opt->avg, opt->std); debug_utig(v, 296); return;
	old_verbose = fm_verbose;
	fm_verbose = 1; // disable all messages and warnings
	t = cputime();
	treal = realtime();
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) {
		w[i].start = i, w[i].step = n_threads;
		w[i].max_dist = max_dist, w[i].opt = opt;
		w[i].e = e, w[i].h = h;
		w[i].v = v;
	}
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], &attr, worker, w + i);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	free(w); free(tid);
	fm_verbose = old_verbose;
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] patched gaps in %.3f sec (%.3f wall-clock sec)\n", __func__, cputime() - t, realtime() - treal);

	if (opt->pr_links)
		for (i = 0; i < v->n; ++i) {
			debug_utig(v, i<<1|0);
			debug_utig(v, i<<1|1);
		}
	make_scaftigs(v, opt->a_thres, opt->p_thres);
	kh_destroy(64, h);
	utig_destroy(v);
}
