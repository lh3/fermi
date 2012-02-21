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
	uint64_t k[2];
	int len, nsr, maxo;
	uint8_t *seq;
	ku128_v reads;
	uint64_t dist[2][2];
	int64_t nei[2];
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
		p->len = kseq->seq.l;
		p->seq = (uint8_t*)strdup(kseq->seq.s);
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

static void collect_nei(utig_v *v, double avg, double std)
{
	int i, j, a, is_absent, max_dist;
	hash64_t *h, *t;
	khint_t k;

	max_dist = (int)(avg + 2. * std + .499);
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
				if (kh_val(t, k) >= p->dist[a][0])
					p->dist[a][1] = p->dist[a][0], p->dist[a][0] = kh_val(t, k), p->nei[a] = kh_key(t, k);
				else if (kh_val(t, k) >= p->dist[a][1]) p->dist[a][1] = kh_val(t, k);
			}
		}
	}
	kh_destroy(64, t);
	kh_destroy(64, h);
#if 1
	for (i = 0; i < v->n; ++i) {
		utig_t *p = &v->a[i];
		for (a = 0; a < 2; ++a)
			if (p->nei[a] >= 0)
				fprintf(stderr, "%d[%ld:%ld]\t%ld\t%d:%ld\t%d:%ld\n", i<<1|a, (long)p->k[0], (long)p->k[1], (long)p->nei[a],
						(int)(p->dist[a][0]>>40), (long)(p->dist[a][0]<<24>>24), (int)(p->dist[a][1]>>40), (long)(p->dist[a][1]<<24>>24));
	}
#endif
}

void mag_scaf_core(const char *fn, double avg, double std, int min_supp)
{
	utig_v *v;
	double rdist;
	v = read_utig(fn, min_supp);
	collect_nei(v, avg, std);
	rdist = cal_rdist(v);
	printf("%f\n", rdist);
}
