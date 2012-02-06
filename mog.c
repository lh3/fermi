/* remaining problems:

  1. multiedges due to tandem repeats
*/

#include <math.h>
#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include "mog.h"
#include "priv.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_INIT2(64,, khint64_t, uint64_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

typedef khash_t(64) hash64_t;

#define ku128_xlt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y > (b).y))
#define ku128_ylt(a, b) ((int64_t)(a).y > (int64_t)(b).y)
#include "ksort.h"
KSORT_INIT(128x, ku128_t, ku128_xlt)
KSORT_INIT(128y, ku128_t, ku128_ylt)
KSORT_INIT_GENERIC(uint64_t)

#define edge_mark_del(_x) ((_x).x = (uint64_t)-2, (_x).y = 0)
#define edge_is_del(_x)   ((_x).x == (uint64_t)-2 || (_x).y == 0)

/*********************
 * Vector operations *
 *********************/

static inline void v128_clean(ku128_v *r)
{
	int i, j;
	for (i = j = 0; i < r->n; ++i)
		if (!edge_is_del(r->a[i])) { // keep this arc
			if (j != i) r->a[j++] = r->a[i];
			else ++j;
		}
	r->n = j;
}

static inline void v128_rmdup(ku128_v *r)
{
	int l, cnt;
	uint64_t x;
	if (r->n > 1) ks_introsort(128x, r->n, r->a);
	for (l = cnt = 0; l < r->n; ++l) // jump to the first node to be retained
		if (edge_is_del(r->a[l])) ++cnt;
		else break;
	if (l == r->n) { // no good arcs
		r->n = 0;
		return;
	}
	x = r->a[l].x;
	for (++l; l < r->n; ++l) { // mark duplicated node
		if (edge_is_del(r->a[l]) || r->a[l].x == x)
			edge_mark_del(r->a[l]), ++cnt;
		else x = r->a[l].x;
	}
	if (cnt) v128_clean(r);
}

static inline void v128_cap(ku128_v *r, int max)
{
	int i, thres;
	if (r->n < max) return;
	ks_introsort(128y, r->n, r->a);
	thres = r->a[max].y;
	for (i = 0; i < r->n; ++i)
		if (r->a[i].y == thres) break;
	r->n = i;
}

/*************************************************
 * Mapping between vertex id and interval end id *
 *************************************************/

static hash64_t *build_hash(const mogv_v *nodes)
{
	long i;
	int j, ret;
	hash64_t *h;
	h = kh_init(64);
	for (i = 0; i < nodes->n; ++i) {
		const mogv_t *p = &nodes->a[i];
		for (j = 0; j < 2; ++j) {
			khint_t k = kh_put(64, h, p->k[j], &ret);
			if (ret == 0) {
				if (fm_verbose >= 2)
					fprintf(stderr, "[W::%s] terminal %ld is duplicated.\n", __func__, (long)p->k[j]);
				kh_val(h, k) = (uint64_t)-1;
			} else kh_val(h, k) = i<<1|j;
		}
	}
	return h;
}

static inline uint64_t tid2idd(hash64_t *h, uint64_t tid)
{
	khint_t k = kh_get(64, h, tid);
	assert(k != kh_end(h));
	return k == kh_end(h)? (uint64_t)-1 : kh_val(h, k);
}

void mog_amend(mog_t *g)
{
	int i, j, l, ll;
	for (i = 0; i < g->v.n; ++i) {
		mogv_t *p = &g->v.a[i];
		ku128_v *r;
		for (j = 0; j < 2; ++j) {
			for (l = 0; l < p->nei[j].n; ++l) {
				uint64_t x = p->nei[j].a[l].x;
				uint64_t z = tid2idd(g->h, x);
				if (z == (uint64_t)-1) { // neighbor is not in the hash table; likely due to tip removal
					edge_mark_del(p->nei[j].a[l]);
					continue;
				}
				r = &g->v.a[z>>1].nei[z&1];
				for (ll = 0; ll < r->n; ++ll)
					if (r->a[ll].x == p->k[j]) break;
				if (ll == r->n) { // not in neighbor's neighor
					p->nei[j].a[l].x = (uint64_t)-1;
					continue;
				}
			}
			v128_rmdup(&p->nei[j]);
		}
	}
}

/*********************************
 * Graph I/O initialization etc. *
 *********************************/

mogopt_t *mog_init_opt()
{
	mogopt_t *o;
	o = calloc(1, sizeof(mogopt_t));
	o->flag |= MOG_F_DROP_TIP0;
	o->max_arc = 512;
	o->min_el = 300;
	o->min_dratio0 = 0.7;
	return o;
}

void mog_v_write(const mogv_t *p, kstring_t *out)
{
	int j, k;
	if (p->len <= 0) return;
	out->l = 0;
	kputc('@', out); kputl(p->k[0], out); kputc(':', out); kputl(p->k[1], out);
	kputc('\t', out); kputw(p->nsr, out);
	for (j = 0; j < 2; ++j) {
		const ku128_v *r = &p->nei[j];
		kputc('\t', out);
		for (k = 0; k < r->n; ++k) {
			if (edge_is_del(r->a[k])) continue;
			kputl(r->a[k].x, out); kputc(',', out); kputw((int32_t)r->a[k].y, out);
			kputc(';', out);
		}
		if (p->nei[j].n == 0) kputc('.', out);
	}
	kputc('\n', out);
	ks_resize(out, out->l + 2 * p->len + 5);
	for (j = 0; j < p->len; ++j)
		out->s[out->l++] = "ACGT"[(int)p->seq[j] - 1];
	out->s[out->l] = 0;
	kputsn("\n+\n", 3, out);
	kputsn(p->cov, p->len, out);
	kputc('\n', out);
}

void mog_g_print(const mog_t *g)
{
	int i;
	kstring_t out;
	out.l = out.m = 0; out.s = 0;
	for (i = 0; i < g->v.n; ++i) {
		if (g->v.a[i].len < 0) continue;
		mog_v_write(&g->v.a[i], &out);
		fwrite(out.s, 1, out.l, stdout);
	}
}

mog_t *mog_g_read(const char *fn, const mogopt_t *opt)
{
	gzFile fp;
	kseq_t *seq;
	ku128_v nei;
	mog_t *g;

	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	kv_init(nei);
	g = calloc(1, sizeof(mog_t));
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i, j;
		char *q;
		mogv_t *p;
		kv_pushp(mogv_t, g->v, &p);
		memset(p, 0, sizeof(mogv_t));
		p->len = -1;
		// parse ->k[2]
		p->k[0] = strtol(seq->name.s, &q, 10); ++q;
		p->k[1] = strtol(q, &q, 10);
		// parse ->nsr
		p->nsr = strtol(seq->comment.s, &q, 10); ++q;
		// parse ->nei[2]
		for (j = 0; j < 2; ++j) {
			int max, max2;
			double thres; // threshold for dropping an overlap
			max = max2 = 0; // largest and 2nd largest overlaps
			nei.n = 0;
			if (*q == '.') {
				q += 2; // skip "." and "\t" (and perhaps "\0", but does not matter)
				continue;
			}
			while (isdigit(*q)) { // parse the neighbors
				ku128_t *r;
				kv_pushp(ku128_t, nei, &r);
				r->x = strtol(q, &q, 10); ++q;
				r->y = strtol(q, &q, 10); ++q;
				g->min_ovlp = g->min_ovlp < r->y? g->min_ovlp : r->y;
				if (max < r->y) max = max2, max = r->y;
				else if (max2 < r->y) max2 = r->y;
			}
			++q; // skip the tailing blank
			thres = (int)(max2 * opt->min_dratio0 + .499);
			for (i = 0; i < nei.n; ++i)
				if (nei.a[i].y < thres) nei.a[i].y = 0; // to be deleted in rmdup_128v()
			v128_rmdup(&nei);
			v128_cap(&nei, opt->max_arc);
			kv_copy(ku128_t, p->nei[j], nei);
		}
		// test if to cut a tip
		p->len = seq->seq.l;
		if (opt->flag & MOG_F_DROP_TIP0) {
			if ((p->nei[0].n == 0 || p->nei[1].n == 0) && p->len < opt->min_el && p->nsr == 1) {
				free(p->nei[0].a); free(p->nei[1].a); // only ->nei[2] have been allocated so far
				--g->v.n;
				continue;
			}
		}
		// set ->{seq,cov,max_len}
		p->max_len = p->len + 1;
		kroundup32(p->max_len);
		p->seq = malloc(p->max_len);
		for (i = 0; i < p->len; ++i) p->seq[i] = seq_nt6_table[(int)seq->seq.s[i]];
		p->cov = malloc(p->max_len);
		strcpy(p->cov, seq->qual.s);
		p->aux[0] = p->aux[1] = -1;
	}
	// free and finalize the graph
	kseq_destroy(seq);
	gzclose(fp);
	free(nei.a);
	g->h = build_hash(&g->v);
	mog_amend(g);
	//g->rdist = fmg_compute_rdist(&g->v);
	return g;
}

/**************************
 * Basic graph operations *
 **************************/

void mog_v_destroy(mogv_t *v)
{
	free(v->nei[0].a); free(v->nei[1].a);
	free(v->seq); free(v->cov);
	memset(v, 0, sizeof(mogv_t));
	v->len = -1;
}

void mog_v_copy_to_empty(mogv_t *dst, const mogv_t *src) // NB: memory leak if dst is allocated
{
	memcpy(dst, src, sizeof(mogv_t));
	dst->max_len = dst->len + 1;
	kroundup32(dst->max_len);
	dst->seq = calloc(dst->max_len, 1); memcpy(dst->seq, src->seq, src->len);
	dst->cov = calloc(dst->max_len, 1); memcpy(dst->cov, src->cov, src->len);
	kv_init(dst->nei[0]); kv_copy(ku128_t, dst->nei[0], src->nei[0]);
	kv_init(dst->nei[1]); kv_copy(ku128_t, dst->nei[1], src->nei[1]);
}

void mog_eh_add(mog_t *g, uint64_t u, uint64_t v, int ovlp) // add v to u
{
	ku128_v *r;
	ku128_t *q;
	uint64_t idd = tid2idd(g->h, u);
	r = &g->v.a[idd>>1].nei[idd&1];
	kv_pushp(ku128_t, *r, &q);
	q->x = v; q->y = ovlp;
}

void mog_eh_markdel(mog_t *g, uint64_t u, uint64_t v) // mark deletion of v from u
{
	int i;	
	uint64_t idd = tid2idd(g->h, u);
	ku128_v *r = &g->v.a[idd>>1].nei[idd&1];
	for (i = 0; i < r->n; ++i)
		if (r->a[i].x == v) edge_mark_del(r->a[i]);
}

void mog_v_del(mog_t *g, mogv_t *p)
{
	int i, j;
	khint_t k;
	for (i = 0; i < 2; ++i) {
		ku128_v *r = &p->nei[i];
		for (j = 0; j < r->n; ++j)
			if (!edge_is_del(r->a[j]))
				mog_eh_markdel(g, r->a[j].x, p->k[i]);
		k = kh_get(64, g->h, p->k[i]);
		kh_del(64, g->h, k);
	}
	mog_v_destroy(p);
}

void mog_v_transdel(mog_t *g, mogv_t *p, int min_ovlp)
{
	if (p->nei[0].n && p->nei[1].n) {
		int i, j, ovlp;
		for (i = 0; i < p->nei[0].n; ++i) {
			if (edge_is_del(p->nei[0].a[i]) || p->nei[0].a[i].x == p->k[0] || p->nei[0].a[i].x == p->k[1]) continue; // due to p->p loop
			for (j = 0; j < p->nei[1].n; ++j) {
				if (edge_is_del(p->nei[1].a[i]) || p->nei[1].a[j].x == p->k[0] || p->nei[1].a[j].x == p->k[1]) continue;
				ovlp = (int)(p->nei[0].a[i].y + p->nei[1].a[j].y) - p->len;
				if (ovlp >= min_ovlp) {
					mog_eh_add(g, p->nei[0].a[i].x, p->nei[1].a[j].x, ovlp);
					mog_eh_add(g, p->nei[1].a[j].x, p->nei[0].a[i].x, ovlp);
				}
			}
		}
	}
	mog_v_del(g, p);
}

#define __swap(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))

void mog_v_flip(mog_t *g, mogv_t *p)
{
	ku128_v t;
	khint_t k;
	hash64_t *h = (hash64_t*)g->h;

	seq_revcomp6(p->len, (uint8_t*)p->seq);
	seq_reverse(p->len, (uint8_t*)p->cov);
	__swap(p->k[0], p->k[1]);
	t = p->nei[0]; p->nei[0] = p->nei[1]; p->nei[1] = t;
	k = kh_get(64, g->h, p->k[0]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
	k = kh_get(64, h, p->k[1]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
}

/*********************
 * Unambiguous merge *
 *********************/

int mog_vh_merge_try(mog_t *g, mogv_t *p) // merge p's neighbor to the right-end of p
{
	mogv_t *q;
	khint_t kp, kq;
	int i, j, new_l;
	hash64_t *h = (hash64_t*)g->h;

	// check if an unambiguous merge can be performed
	if (p->nei[1].n != 1) return -1; // multiple or no neighbor; do not merge
	kq = kh_get(64, g->h, p->nei[1].a[0].x); assert(kq != kh_end(h)); // otherwise the neighbor is non-existant
	q = &g->v.a[kh_val((hash64_t*)g->h, kq)>>1];
	if (p == q) return -2; // we have a loop p->p. We cannot merge in this case
	if (q->nei[kh_val(h, kq)&1].n != 1) return -3; // the neighbor q has multiple neighbors. cannot be an unambiguous merge

	// we can perform a merge; do further consistency check (mostly check bugs)
	if (kh_val(h, kq)&1) mog_v_flip(g, q); // a "><" bidirectional arc; flip q
	kp = kh_get(64, g->h, p->k[1]); assert(kp != kh_end(h)); // get the iterator to p
	kh_del(64, g->h, kp); kh_del(64, g->h, kq); // remove the two ends of the arc in the hash table
	assert(p->k[1] == q->nei[0].a[0].x && q->k[0] == p->nei[1].a[0].x); // otherwise inconsistent topology
	assert(p->nei[1].a[0].y == q->nei[0].a[0].y); // the overlap length must be the same
	assert(p->len >= p->nei[1].a[0].y && q->len >= p->nei[1].a[0].y); // and the overlap is shorter than both vertices

	// update the read count and sequence length
	p->nsr += q->nsr;
	new_l = p->len + q->len - p->nei[1].a[0].y;
	if (new_l + 1 > p->max_len) { // then double p->seq and p->cov
		p->max_len = new_l + 1;
		kroundup32(p->max_len);
		p->seq = realloc(p->seq, p->max_len);
		p->cov = realloc(p->cov, p->max_len);
	}
	// merge seq and cov
	for (i = p->len - p->nei[1].a[0].y, j = 0; j < q->len; ++i, ++j) { // write seq and cov
		p->seq[i] = q->seq[j];
		if (i < p->len) {
			if ((int)p->cov[i] + (q->cov[j] - 33) > 126) p->cov[i] = 126;
			else p->cov[i] += q->cov[j] - 33;
		} else p->cov[i] = q->cov[j];
	}
	p->seq[new_l] = p->cov[new_l] = 0;
	p->len = new_l;
	// merge neighbors
	free(p->nei[1].a);
	p->nei[1] = q->nei[1]; p->k[1] = q->k[1];
	q->nei[1].a = 0; // to avoid freeing p->nei[1] by mog_v_destroy() below
	// update the hash table for the right end of p
	kp = kh_get(64, g->h, p->k[1]);
	assert(kp != kh_end((hash64_t*)g->h));
	kh_val(h, kp) = (p - g->v.a)<<1 | 1;
	// clean up q
	mog_v_destroy(q);
	return 0;
}

void mog_g_merge(mog_t *g)
{
	int i;
	double tcpu = cputime();
	for (i = 0; i < g->v.n; ++i) { // remove multiedges; FIXME: should we do that?
		v128_rmdup(&g->v.a[i].nei[0]);
		v128_rmdup(&g->v.a[i].nei[1]);
	}
	for (i = 0; i < g->v.n; ++i) {
		mogv_t *p = &g->v.a[i];
		if (p->len < 0) continue;
		while (mog_vh_merge_try(g, p) == 0);
		mog_v_flip(g, p);
		while (mog_vh_merge_try(g, p) == 0);
	}
	if (fm_verbose >= 2)
		fprintf(stderr, "[M::%s] merged unambiguous arcs in %.2f sec\n", __func__, cputime() - tcpu);
}

/*****************************
 * Easy graph simplification *
 *****************************/

void mog_g_rm_vext(mog_t *g, int min_len, int min_nsr)
{
	int i;
	for (i = 0; i < g->v.n; ++i) {
		mogv_t *p = &g->v.a[i];
		if ((p->nei[0].n == 0 || p->nei[1].n == 0) && p->len < min_len && p->nsr < min_nsr)
			mog_v_del(g, p);
	}
}

void mog_g_rm_vint(mog_t *g, int min_len, int min_nsr, int min_ovlp)
{
	int i;
	for (i = 0; i < g->v.n; ++i) {
		mogv_t *p = &g->v.a[i];
		if (p->len < min_len && p->nsr < min_nsr)
			mog_v_transdel(g, p, min_ovlp);
	}
}

void mog_g_rm_edge(mog_t *g, int min_ovlp, double min_ratio)
{
	int i, j, k;
	for (i = 0; i < g->v.n; ++i) {
		mogv_t *p = &g->v.a[i];
		for (j = 0; j < 2; ++j) {
			ku128_v *r = &p->nei[j];
			int max_ovlp = min_ovlp;
			if (r->n == 0) continue; // no overlapping reads
			for (k = 0; k < r->n; ++k) // get the max overlap length
				if (max_ovlp < r->a[k].y) max_ovlp = r->a[k].y;
			for (k = 0; k < r->n; ++k) {
				if (r->a[k].y < min_ovlp || (double)r->a[k].y / max_ovlp < min_ratio) {
					mog_eh_markdel(g, r->a[k].x, p->k[j]); // FIXME: should we check if r->a[k] is p itself?
					edge_mark_del(r->a[k]);
				}
			}
		}
	}
}

/*********************************************
 * A-statistics and simplistic flow analysis *
 *********************************************/

#define A_THRES 20.
#define A_MIN_SUPP 5

double mog_cal_rdist(mog_t *g)
{
	mogv_v *v = &g->v;
	int j;
	uint64_t *srt;
	double rdist = -1., t;
	int64_t i, sum_n_all, sum_n, sum_l;

	t = cputime();
	srt = calloc(v->n, 8);
	for (i = 0, sum_n_all = 0; i < v->n; ++i) {
		srt[i] = (uint64_t)v->a[i].nsr<<32 | i;
		sum_n_all += v->a[i].nsr;
	}
	ks_introsort_uint64_t(v->n, srt);

	for (j = 0; j < 2; ++j) {
		sum_n = sum_l = 0;
		for (i = v->n - 1; i >= 0; --i) {
			const mogv_t *p = &v->a[srt[i]<<32>>32];
			int tmp1, tmp2;
			tmp1 = tmp2 = 0;
			if (p->nei[0].n) ++tmp1, tmp2 += p->nei[0].a[0].y;
			if (p->nei[1].n) ++tmp1, tmp2 += p->nei[1].a[0].y;
			if (tmp1) tmp2 /= tmp1;
			if (rdist > 0.) {
				double A = (p->len - tmp1) / rdist - p->nsr * M_LN2;
				if (A < A_THRES) continue;
			}
			sum_n += p->nsr;
			sum_l += p->len - tmp1;
			if (sum_n >= sum_n_all * 0.5) break;
		}
		rdist = (double)sum_l / sum_n;
	}
	if (fm_verbose >= 3) {
		fprintf(stderr, "[M::%s] average read distance %.3f, computed in %.3f seconds.\n", __func__, rdist, cputime() - t);
		fprintf(stderr, "[M::%s] approximate genome size: %.0f (inaccurate!)\n", __func__, rdist * sum_n_all);
	}

	free(srt);
	return rdist;
}

void mog_vh_flowflt(mog_t *g, size_t idd, double thres)
{ // only works for: {p,r}->q, where both p and q are unique and p has a single neighbor q
	mogv_t *p, *q, *t;
	double A;
	uint64_t u, v;
	ku128_v *r;
	int i;

	p = &g->v.a[idd>>1];
	if (p->nei[idd&1].n != 1) return;
	A = (p->len - p->nei[idd&1].a[0].y) / g->rdist - p->nsr * M_LN2;
	if (A < thres && p->nsr < A_MIN_SUPP) return; // p not significantly unique
	u = tid2idd(g->h, p->nei[idd&1].a[0].x);
	q = &g->v.a[u>>1];
	if (p == q) return;
	if (q->nei[u&1].n < 2) return; // well, p and q can be merged already
	assert(q->len >= p->nei[idd&1].a[0].y);
	A = (q->len - p->nei[idd&1].a[0].y) / g->rdist - q->nsr * M_LN2;
	if (A < thres) return; // q not significantly unique
//	fprintf(stderr, "%lld:%lld, %lld:%lld\n", p->k[0], p->k[1], q->k[0], q->k[1]);

	r = &q->nei[u&1];
	for (i = 0; i < r->n; ++i) {
		int to_cut = 0;
		v = tid2idd(g->h, r->a[i].x);
		if (v == (uint64_t)-1) continue;
		t = &g->v.a[v>>1];
		if (t->nei[v&1].n == 1) {
			// Should we also consider to cut in this case? Perhaps does not matter.
			if (t->nsr <= 2 && t->nei[(v&1)^1].n == 0) to_cut = 1; // a small tip
		} else to_cut = 1;
		if (to_cut) {
			if (r->a[i].x != q->k[0] && r->a[i].x != q->k[1])
				mog_eh_markdel(g, r->a[i].x, q->k[u&1]);
			edge_mark_del(r->a[i]);
		}
	}
}

/*************************
 * SW high-coverage tips *
 *************************/

#include "ksw.h"

void mog_v_swrm(mog_t *g, mogv_t *p)
{
	int i, j, k, l, dir, max_l, l_qry;
	mogv_t *q, *t;
	ku128_v *r, *s;
	uint8_t *seq;
	int8_t mat[16];
	ksw_query_t *qry;
	ksw_aux_t aux;

	if (p->len < 0) return;
	//if (p->nei[0].n && p->nei[1].n) return; // FIXME: between this and the next line, which is better?
	if (p->nei[0].n + p->nei[1].n != 1) return;
	dir = p->nei[0].n? 0 : 1;
	// initialize the scoring system
	for (i = k = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? 5 : -4;
	aux.gapo = 6; aux.gape = 3;
	
	s = &p->nei[dir];
	for (l = 0; l < s->n; ++l) { // if we use "if (p->nei[0].n + p->nei[1].n != 1)", s->n == 1
		uint64_t v;
		v = tid2idd(g->h, s->a[l].x);
		q = &g->v.a[v>>1];
		if (q == p || q->nei[v&1].n == 1) continue;
		// get the query ready
		max_l = (p->len - s->a[l].y) * 2;
		seq = malloc(max_l + 1);
		if (dir == 0) { // forward strand
			for (j = s->a[l].y, k = 0; j < p->len; ++j)
				seq[k++] = p->seq[j] - 1;
		} else { // reverse
			for (j = p->len - s->a[l].y - 1, k = 0; j >= 0; --j)
				seq[k++] = 4 - p->seq[j];
		}
		l_qry = k; aux.T = l_qry / 2;
		qry = ksw_qinit(2, l_qry, seq, 4, mat);
		//fprintf(stderr, "===> %lld:%lld:%d[%d], %d, %ld <===\n", p->k[0], p->k[1], s->n, l, p->n, q->nei[v&1].n);
		//for (j = 0; j < k; ++j) fputc("ACGTN"[(int)seq[j]], stderr); fputc('\n', stderr);

		r = &q->nei[v&1];
		for (i = 0; i < r->n; ++i) {
			uint64_t w;
			if (r->a[i].x == p->k[dir]) continue;
			w = tid2idd(g->h, r->a[i].x);
			// get the target sequence
			t = &g->v.a[w>>1];
			if (w&1) { // reverse strand
				for (j = t->len - r->a[i].y - 1, k = 0; j >= 0 && k < max_l; --j)
					seq[k++] = 4 - t->seq[j];
			} else {
				for (j = r->a[i].y, k = 0; j < t->len && k < max_l; ++j)
					seq[k++] = t->seq[j] - 1;
			}
			ksw_sse2(qry, k, seq, &aux);
			//for (j = 0; j < k; ++j) fputc("ACGTN"[(int)seq[j]], stderr); fprintf(stderr, "\t%d\n", aux.score);
			if (aux.score) {
				double r_diff, n_diff;
				n_diff = (l_qry * 5. - aux.score) / (5. + 4.); // 5: matching score; -4: mismatchig score
				r_diff = n_diff / l_qry;
				if ((int)(n_diff + .499) <= 1 || r_diff < 0.1) break;
			}
		}

		if (i != r->n) {
			// mark delete in p and delete in q
			edge_mark_del(s->a[l]);
			for (i = 0; i < r->n; ++i)
				if (r->a[i].x != p->k[dir])
					edge_mark_del(r->a[i]);
		}
		free(seq); free(qry);
	}

	for (i = 0; i < s->n; ++i)
		if (!edge_is_del(s->a[i])) break;
	if (i == s->n) mog_v_del(g, p); // p is not connected to any other vertices
}
