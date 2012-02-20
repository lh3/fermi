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

#define A_THRES 20.

extern unsigned char seq_nt6_table[128];

typedef struct {
	uint64_t k[2];
	int len, nsr, maxo;
	uint8_t *seq;
	ku128_v reads;
} utig_t;

typedef kvec_t(utig_t) utig_v;

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

		kv_init(p->reads);
		while (isdigit(*q)) { // read mapping
			ku128_t x;
			x.x = strtol(q, &q, 10); ++q;
			x.y = x.x&1; x.x >>= 1; // to fit the msg format; FIXME: better make this cleaner...
			x.y |= (uint64_t)strtol(q, &q, 10)<<32; ++q;
			x.y |= strtol(q, &q, 10)<<1;
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

void mag_scaf_core(const char *fn, double avg, double std, int min_supp)
{
	utig_v *v;
	double rdist;
	v = read_utig(fn, min_supp);
	rdist = cal_rdist(v);
	printf("%f\n", rdist);
}
