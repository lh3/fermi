/* remaining problems:

  1. multiedges due to tandem repeats
*/

#include <math.h>
#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include "mag.h"
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

void mag_v128_clean(ku128_v *r)
{
	v128_clean(r);
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
	if (r->n <= max) return;
	ks_introsort(128y, r->n, r->a);
	thres = r->a[max].y;
	for (i = 0; i < r->n; ++i)
		if (r->a[i].y == thres) break;
	r->n = i;
}

/*************************************************
 * Mapping between vertex id and interval end id *
 *************************************************/

void mag_g_build_hash(mag_t *g)
{
	long i;
	int j, ret;
	hash64_t *h;
	h = kh_init(64);
	for (i = 0; i < g->v.n; ++i) {
		const magv_t *p = &g->v.a[i];
		for (j = 0; j < 2; ++j) {
			khint_t k = kh_put(64, h, p->k[j], &ret);
			if (ret == 0) {
				if (fm_verbose >= 2)
					fprintf(stderr, "[W::%s] terminal %ld is duplicated.\n", __func__, (long)p->k[j]);
				kh_val(h, k) = (uint64_t)-1;
			} else kh_val(h, k) = i<<1|j;
		}
	}
	g->h = h;
}

static inline uint64_t tid2idd(hash64_t *h, uint64_t tid)
{
	khint_t k = kh_get(64, h, tid);
	assert(k != kh_end(h));
	return kh_val(h, k);
}

uint64_t mag_tid2idd(void *h, uint64_t tid) // exported version
{
	return tid2idd(h, tid);
}

void mag_amend(mag_t *g)
{
	int i, j, l, ll;
	for (i = 0; i < g->v.n; ++i) {
		magv_t *p = &g->v.a[i];
		ku128_v *r;
		for (j = 0; j < 2; ++j) {
			for (l = 0; l < p->nei[j].n; ++l) {
				khint_t k;
				uint64_t z, x = p->nei[j].a[l].x;
				k = kh_get(64, g->h, x);
				if (k == kh_end((hash64_t*)g->h)) { // neighbor is not in the hash table; likely due to tip removal
					edge_mark_del(p->nei[j].a[l]);
					continue;
				} else z = kh_val((hash64_t*)g->h, k);
				r = &g->v.a[z>>1].nei[z&1];
				for (ll = 0, z = p->k[j]; ll < r->n; ++ll)
					if (r->a[ll].x == z) break;
				if (ll == r->n) // not in neighbor's neighor
					edge_mark_del(p->nei[j].a[l]);
			}
			v128_rmdup(&p->nei[j]);
		}
	}
}

/*********************************
 * Graph I/O initialization etc. *
 *********************************/

void mag_v_write(const magv_t *p, kstring_t *out)
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

void mag_g_print(const mag_t *g)
{
	int i;
	kstring_t out;
	out.l = out.m = 0; out.s = 0;
	for (i = 0; i < g->v.n; ++i) {
		if (g->v.a[i].len < 0) continue;
		mag_v_write(&g->v.a[i], &out);
		fwrite(out.s, 1, out.l, stdout);
	}
	free(out.s);
	fflush(stdout);
}

mag_t *mag_g_read(const char *fn, const magopt_t *opt)
{
	gzFile fp;
	kseq_t *seq;
	ku128_v nei;
	mag_t *g;
	int is_mod = 0;
	double t;

	t = cputime();
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	kv_init(nei);
	g = calloc(1, sizeof(mag_t));
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i, j;
		char *q;
		magv_t *p;
		kv_pushp(magv_t, g->v, &p);
		memset(p, 0, sizeof(magv_t));
		p->len = -1;
		// parse ->k[2]
		p->k[0] = strtol(seq->name.s, &q, 10); ++q;
		p->k[1] = strtol(q, &q, 10);
		// parse ->nsr
		p->nsr = strtol(seq->comment.s, &q, 10); ++q;
		// parse ->nei[2]
		for (j = 0; j < 2; ++j) {
			int max, max2;
			max = max2 = 0; // largest and 2nd largest overlaps
			nei.n = 0;
			if (*q == '.') {
				q += 2; // skip "." and "\t" (and perhaps "\0", but does not matter)
				continue;
			}
			while (isdigit(*q) || *q == '-') { // parse the neighbors
				ku128_t *r;
				kv_pushp(ku128_t, nei, &r);
				r->x = strtol(q, &q, 10); ++q;
				r->y = strtol(q, &q, 10); ++q;
				g->min_ovlp = g->min_ovlp < r->y? g->min_ovlp : r->y;
				if (max < r->y) max = max2, max = r->y;
				else if (max2 < r->y) max2 = r->y;
			}
			++q; // skip the tailing blank
			if (!(opt->flag & MOG_F_READ_ORI)) {
				double thres = (int)(max2 * opt->min_dratio0 + .499);
				for (i = 0; i < nei.n; ++i)
					if (nei.a[i].y < thres) is_mod = 1, nei.a[i].y = 0; // to be deleted in rmdup_128v()
				v128_rmdup(&nei);
				if (nei.n > opt->max_arc) {
					is_mod = 1;
					v128_cap(&nei, opt->max_arc);
				}
			}
			kv_copy(ku128_t, p->nei[j], nei);
		}
		// test if to cut a tip
		p->len = seq->seq.l;
		if (!(opt->flag & MOG_F_READ_ORI) && (p->nei[0].n == 0 || p->nei[1].n == 0) && p->len < opt->min_elen && p->nsr == 1) {
			free(p->nei[0].a); free(p->nei[1].a); // only ->nei[2] have been allocated so far
			--g->v.n;
			is_mod = 1;
			continue;
		}
		// set ->{seq,cov,max_len}
		p->max_len = p->len + 1;
		kroundup32(p->max_len);
		p->seq = malloc(p->max_len);
		for (i = 0; i < p->len; ++i) p->seq[i] = seq_nt6_table[(int)seq->seq.s[i]];
		p->cov = malloc(p->max_len);
		if (seq->qual.l == 0)
			for (i = 0; i < seq->seq.l; ++i)
				p->cov[i] = 34;
		else strcpy(p->cov, seq->qual.s);
	}
	// free and finalize the graph
	kseq_destroy(seq);
	gzclose(fp);
	free(nei.a);
	// finalize
	mag_g_build_hash(g);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] read the graph and constructed the dictionary in %.3f sec\n", __func__, cputime() - t);
	if (is_mod && fm_verbose >= 3)
		fprintf(stderr, "[M::%s] the graph is modified during reading.\n", __func__);
	if (is_mod || !(opt->flag & MOG_F_NO_AMEND)) {
		t = cputime();
		mag_amend(g);
		fprintf(stderr, "[M::%s] amended the graph in %.3f sec.\n", __func__, cputime() - t);
	}
	g->rdist = mag_cal_rdist(g);
	if (opt->flag & MOG_F_READnMERGE) mag_g_merge(g, 1);
	return g;
}

/**************************
 * Basic graph operations *
 **************************/

void mag_v_destroy(magv_t *v)
{
	free(v->nei[0].a); free(v->nei[1].a);
	free(v->seq); free(v->cov);
	memset(v, 0, sizeof(magv_t));
	v->len = -1;
}

void mag_g_destroy(mag_t *g)
{
	int i;
	kh_destroy(64, g->h);
	for (i = 0; i < g->v.n; ++i)
		mag_v_destroy(&g->v.a[i]);
	free(g->v.a);
	free(g);
}

void mag_v_copy_to_empty(magv_t *dst, const magv_t *src) // NB: memory leak if dst is allocated
{
	memcpy(dst, src, sizeof(magv_t));
	dst->max_len = dst->len + 1;
	kroundup32(dst->max_len);
	dst->seq = calloc(dst->max_len, 1); memcpy(dst->seq, src->seq, src->len);
	dst->cov = calloc(dst->max_len, 1); memcpy(dst->cov, src->cov, src->len);
	kv_init(dst->nei[0]); kv_copy(ku128_t, dst->nei[0], src->nei[0]);
	kv_init(dst->nei[1]); kv_copy(ku128_t, dst->nei[1], src->nei[1]);
}

void mag_eh_add(mag_t *g, uint64_t u, uint64_t v, int ovlp) // add v to u
{
	ku128_v *r;
	ku128_t *q;
	uint64_t idd;
	int i;
	if ((int64_t)u < 0) return;
	idd = tid2idd(g->h, u);
	r = &g->v.a[idd>>1].nei[idd&1];
	for (i = 0; i < r->n; ++i) // no multi-edges
		if (r->a[i].x == v) return;
	kv_pushp(ku128_t, *r, &q);
	q->x = v; q->y = ovlp;
}

void mag_eh_markdel(mag_t *g, uint64_t u, uint64_t v) // mark deletion of v from u
{
	int i;	
	uint64_t idd;
	if ((int64_t)u < 0) return;
	idd = tid2idd(g->h, u);
	ku128_v *r = &g->v.a[idd>>1].nei[idd&1];
	for (i = 0; i < r->n; ++i)
		if (r->a[i].x == v) edge_mark_del(r->a[i]);
}

void mag_v_del(mag_t *g, magv_t *p)
{
	int i, j;
	khint_t k;
	if (p->len < 0) return;
	for (i = 0; i < 2; ++i) {
		ku128_v *r = &p->nei[i];
		for (j = 0; j < r->n; ++j)
			if (!edge_is_del(r->a[j]) && r->a[j].x != p->k[0] && r->a[j].x != p->k[1])
				mag_eh_markdel(g, r->a[j].x, p->k[i]);
	}
	for (i = 0; i < 2; ++i) {
		k = kh_get(64, g->h, p->k[i]);
		kh_del(64, g->h, k);
	}
	mag_v_destroy(p);
}

void mag_v_transdel(mag_t *g, magv_t *p, int min_ovlp)
{
	if (p->nei[0].n && p->nei[1].n) {
		int i, j, ovlp;
		for (i = 0; i < p->nei[0].n; ++i) {
			if (edge_is_del(p->nei[0].a[i]) || p->nei[0].a[i].x == p->k[0] || p->nei[0].a[i].x == p->k[1]) continue; // due to p->p loop
			for (j = 0; j < p->nei[1].n; ++j) {
				if (edge_is_del(p->nei[1].a[j]) || p->nei[1].a[j].x == p->k[0] || p->nei[1].a[j].x == p->k[1]) continue;
				ovlp = (int)(p->nei[0].a[i].y + p->nei[1].a[j].y) - p->len;
				if (ovlp >= min_ovlp) {
					mag_eh_add(g, p->nei[0].a[i].x, p->nei[1].a[j].x, ovlp);
					mag_eh_add(g, p->nei[1].a[j].x, p->nei[0].a[i].x, ovlp);
				}
			}
		}
	}
	mag_v_del(g, p);
}

void mag_v_flip(mag_t *g, magv_t *p)
{
	ku128_v t;
	khint_t k;
	hash64_t *h = (hash64_t*)g->h;

	seq_revcomp6(p->len, (uint8_t*)p->seq);
	seq_reverse(p->len, (uint8_t*)p->cov);
	p->k[0] ^= p->k[1]; p->k[1] ^= p->k[0]; p->k[0] ^= p->k[1];
	t = p->nei[0]; p->nei[0] = p->nei[1]; p->nei[1] = t;
	k = kh_get(64, h, p->k[0]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
	k = kh_get(64, h, p->k[1]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
}

/*********************
 * Unambiguous merge *
 *********************/

int mag_vh_merge_try(mag_t *g, magv_t *p) // merge p's neighbor to the right-end of p
{
	magv_t *q;
	khint_t kp, kq;
	int i, j, new_l;
	hash64_t *h = (hash64_t*)g->h;

	// check if an unambiguous merge can be performed
	if (p->nei[1].n != 1) return -1; // multiple or no neighbor; do not merge
	if ((int64_t)p->nei[1].a[0].x < 0) return -2;
	kq = kh_get(64, g->h, p->nei[1].a[0].x);
	assert(kq != kh_end(h)); // otherwise the neighbor is non-existant
	q = &g->v.a[kh_val((hash64_t*)g->h, kq)>>1];
	if (p == q) return -3; // we have a loop p->p. We cannot merge in this case
	if (q->nei[kh_val(h, kq)&1].n != 1) return -4; // the neighbor q has multiple neighbors. cannot be an unambiguous merge

	// we can perform a merge; do further consistency check (mostly check bugs)
	if (kh_val(h, kq)&1) mag_v_flip(g, q); // a "><" bidirectional arc; flip q
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
	q->nei[1].a = 0; // to avoid freeing p->nei[1] by mag_v_destroy() below
	// update the hash table for the right end of p
	kp = kh_get(64, g->h, p->k[1]);
	assert(kp != kh_end((hash64_t*)g->h));
	kh_val(h, kp) = (p - g->v.a)<<1 | 1;
	// clean up q
	mag_v_destroy(q);
	return 0;
}

void mag_g_merge(mag_t *g, int rmdup)
{
	int i;
	for (i = 0; i < g->v.n; ++i) { // remove multiedges; FIXME: should we do that?
		if (rmdup) {
			v128_rmdup(&g->v.a[i].nei[0]);
			v128_rmdup(&g->v.a[i].nei[1]);
		} else {
			v128_clean(&g->v.a[i].nei[0]);
			v128_clean(&g->v.a[i].nei[1]);
		}
	}
	for (i = 0; i < g->v.n; ++i) {
		magv_t *p = &g->v.a[i];
		if (p->len < 0) continue;
		while (mag_vh_merge_try(g, p) == 0);
		mag_v_flip(g, p);
		while (mag_vh_merge_try(g, p) == 0);
	}
}

/*****************************
 * Easy graph simplification *
 *****************************/

void mag_g_rm_vext(mag_t *g, int min_len, int min_nsr)
{
	int i;
	for (i = 0; i < g->v.n; ++i) {
		magv_t *p = &g->v.a[i];
		if (p->len >= 0 && (p->nei[0].n == 0 || p->nei[1].n == 0) && p->len < min_len && p->nsr < min_nsr)
			mag_v_del(g, p);
	}
}

void mag_g_rm_vint(mag_t *g, int min_len, int min_nsr, int min_ovlp)
{
	int i;
	for (i = 0; i < g->v.n; ++i) {
		magv_t *p = &g->v.a[i];
		if (p->len >= 0 && p->len < min_len && p->nsr < min_nsr)
			mag_v_transdel(g, p, min_ovlp);
	}
}

void mag_g_rm_edge(mag_t *g, int min_ovlp, double min_ratio, int min_len, int min_nsr)
{
	int i, j, k;
	for (i = 0; i < g->v.n; ++i) {
		magv_t *p = &g->v.a[i];
		if (p->len >= 0 && (p->nei[0].n == 0 || p->nei[1].n == 0) && p->len < min_len && p->nsr < min_nsr)
			continue; // skip tips
		for (j = 0; j < 2; ++j) {
			ku128_v *r = &p->nei[j];
			int max_ovlp = min_ovlp, max_k = -1;
			if (r->n == 0) continue; // no overlapping reads
			for (k = 0; k < r->n; ++k) // get the max overlap length
				if (max_ovlp < r->a[k].y)
					max_ovlp = r->a[k].y, max_k = k;
			if (max_k >= 0) { // test if max_k is a tip
				uint64_t x = tid2idd(g->h, r->a[max_k].x);
				magv_t *q = &g->v.a[x>>1];
				if (q->len >= 0 && (q->nei[0].n == 0 || q->nei[1].n == 0) && q->len < min_len && q->nsr < min_nsr)
					max_ovlp = min_ovlp;
			}
			for (k = 0; k < r->n; ++k) {
				if (edge_is_del(r->a[k])) continue;
				if (r->a[k].y < min_ovlp || (double)r->a[k].y / max_ovlp < min_ratio) {
					mag_eh_markdel(g, r->a[k].x, p->k[j]); // FIXME: should we check if r->a[k] is p itself?
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

double mag_cal_rdist(mag_t *g)
{
	magv_v *v = &g->v;
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
			const magv_t *p = &v->a[srt[i]<<32>>32];
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

/**************
 * Key portal *
 **************/

magopt_t *mag_init_opt()
{
	magopt_t *o;

	o = calloc(1, sizeof(magopt_t));
	o->flag = MOG_F_READnMERGE;
	o->max_arc = 512;
	o->min_dratio0 = 0.7;

	o->n_iter = 3;
	o->min_elen = 300;
	o->min_ovlp = 60;
	o->min_ensr = 4;
	o->min_insr = 3;
	o->min_dratio1 = 0.8;

	o->max_bcov = 10.;
	o->max_bfrac = 0.15;
	o->max_bvtx = 64;
	o->max_bdist = 512;
	return o;
}

void mag_g_clean(mag_t *g, const magopt_t *opt)
{
	double t;
	int j;

	if ((opt->flag & MOG_F_CLEAN) == 0) return;
	if (g->min_ovlp < opt->min_ovlp) g->min_ovlp = opt->min_ovlp;
	//mag_vh_simplify_bubble(g, tid2idd(g->h, 34356802), 512, 500, a); exit(0); // a good case
	mag_g_rm_vext(g, opt->min_elen, opt->min_ensr < 3? opt->min_ensr : 3);
	for (j = 0; j < opt->n_iter; ++j) {
		double r = opt->n_iter == 1? 1. : .5 + .5 * j / (opt->n_iter - 1);
		t = cputime();
		mag_g_rm_edge(g, opt->min_ovlp * r, opt->min_dratio1 * r, opt->min_elen, opt->min_ensr);
		mag_g_rm_vext(g, opt->min_elen * r, opt->min_ensr * r > 2.? opt->min_ensr * r > 2. : 2);
		mag_g_merge(g, 1);
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] finished simple graph simplification round %d in %.3f sec.\n", __func__, j+1, cputime() - t);
	}
	t = cputime();
	for (j = 0; j < opt->n_iter; ++j) {
		mag_g_rm_vext(g, opt->min_elen, opt->min_ensr);
		mag_g_merge(g, 0);
	}
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] finished another %d rounds of tip removal in %.3f sec.\n", __func__, opt->n_iter, cputime() - t);
	if (opt->flag & MOG_F_AGGRESSIVE) {
		t = cputime();
		mag_g_pop_open(g, opt->min_elen);
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] popped open bubbles in %.3f sec.\n", __func__, cputime() - t);
	}
	if (!(opt->flag & MOG_F_NO_SIMPL)) {
		t = cputime();
		mag_g_simplify_bubble(g, opt->max_bvtx, opt->max_bdist);
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] simplified complex bubbles in %.3f sec.\n", __func__, cputime() - t);
	}
	t = cputime();
	mag_g_pop_simple(g, opt->max_bcov, opt->max_bfrac, opt->flag & MOG_F_AGGRESSIVE);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] popped closed bubbles in %.3f sec.\n", __func__, cputime() - t);
	if (opt->min_insr >= 2) {
		t = cputime();
		mag_g_rm_vint(g, opt->min_elen, opt->min_insr, g->min_ovlp);
		mag_g_rm_edge(g, opt->min_ovlp, opt->min_dratio1, opt->min_elen, opt->min_ensr);
		mag_g_rm_vext(g, opt->min_elen, opt->min_ensr);
		mag_g_merge(g, 1);
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] removed interval low-cov vertices in %.3f sec.\n", __func__, cputime() - t);
	}
	t = cputime();
	if (opt->flag & MOG_F_AGGRESSIVE) mag_g_pop_open(g, opt->min_elen);
	else {
		mag_g_rm_vext(g, opt->min_elen, opt->min_ensr);
		mag_g_merge(g, 0);
	}
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] coverage based graph cleanup in %.3f sec.\n", __func__, cputime() - t);
}
