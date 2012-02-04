#include <zlib.h>
#include "priv.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_INIT2(64,, khint64_t, uint64_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

#define fm128_xlt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y > (b).y))
#define fm128_ylt(a, b) ((int64_t)(a).y > (int64_t)(b).y)
#include "ksort.h"
KSORT_INIT(128x, fm128_t, fm128_xlt)
KSORT_INIT(128y, fm128_t, fm128_ylt)

#define arc_mark_del(_x) ((_x).x = (uint64_t)-2, (_x).y = 0)
#define arc_is_del(_x)   ((_x).x == (uint64_t)-2 || (_x).y == 0)

/*********************
 * Vector operations *
 *********************/

static inline void v128_clean(mog128_v *r)
{
	int i, j;
	for (i = j = 0; i < r->n; ++i)
		if (!arc_is_del(r->a[i])) { // keep this arc
			if (j != i) r->a[j++] = r->a[i];
			else ++j;
		}
	r->n = j;
}

static inline void v128_rmdup(mog128_v *r)
{
	int l, cnt;
	uint64_t x;
	if (r->n > 1) ks_introsort(128x, r->n, r->a);
	for (l = cnt = 0; l < r->n; ++l) // jump to the first node to be retained
		if (arc_is_del(r->a[l])) ++cnt;
		else break;
	if (l == r->n) { // no good arcs
		r->n = 0;
		return;
	}
	x = r->a[l].x;
	for (++l; l < r->n; ++l) { // mark duplicated node
		if (arc_is_del(r->a[l]) || r->a[l].x == x)
			arc_mark_del(r->a[l]), ++cnt;
		else x = r->a[l].x;
	}
	if (cnt) v128_clean(r);
}

static inline void v128_cap(mog128_v *r, int max)
{
	int i, thres;
	if (r->n < max) return;
	ks_introsort(128y, r->n, r->a);
	thres = r->a[max].y;
	for (i = 0; i < r->n; ++i)
		if (r->a[i].y == thres) break;
	r->n = i;
}

/**************************************
 * Mapping between node id and end id *
 **************************************/

static hash64_t *build_hash(const mognode_v *nodes)
{
	long i;
	int j, ret;
	hash64_t *h;
	h = kh_init(64);
	for (i = 0; i < nodes->n; ++i) {
		const fmnode_t *p = &nodes->a[i];
		for (j = 0; j < 2; ++j) {
			khint_t k = kh_put(64, h, p->k[j], &ret);
			if (ret == 0) {
				if (fm_verbose >= 2)
					fprintf(stderr, "[W::%s] end %ld is duplicated.\n", __func__, (long)p->k[j]);
				kh_val(h, k) = (uint64_t)-1;
			} else kh_val(h, k) = i<<1|j;
		}
	}
	return h;
}

static inline uint64_t tid2idd(hash64_t *h, uint64_t tid)
{
	khint_t k = kh_get(64, h, tid);
	return k == kh_end(h)? (uint64_t)(-1) : kh_val(h, k);
}

void mog_amend(mog_t *g)
{
	size_t i;
	int j, l, ll;
	for (i = 0; i < g->nodes.n; ++i) {
		mognode_t *p = &g->nodes.a[i];
		mog128_v *r;
		for (j = 0; j < 2; ++j) {
			for (l = 0; l < p->nei[j].n; ++l) {
				uint64_t x = p->nei[j].a[l].x;
				uint64_t z = tid2idd(g->h, x);
				if (z == (uint64_t)-1) { // neighbor is not in the hash table; likely due to tip removal
					arc_mark_del(p->nei[j].a[l]);
					continue;
				}
				r = &g->nodes.a[z>>1].nei[z&1];
				for (ll = 0; ll < r->n; ++ll)
					if (r->a[ll].x == p->k[j]) break;
				if (ll == r->n) { // not in neighbor's neighor
					p->nei[j].a[l].x = (uint64_t)-1;
					continue;
				}
			}
			v128_clean(&p->nei[j]);
		}
	}
}

/*************
 * Graph I/O *
 *************/

void mog_write1(const mognode_t *p, kstring_t *out)
{
	int j, k;
	if (p->len <= 0) return;
	out->l = 0;
	kputc('@', out); kputl(p->k[0], out); kputc(':', out); kputl(p->k[1], out);
	kputc('\t', out); kputw(p->nsr, out);
	for (j = 0; j < 2; ++j) {
		kputc('\t', out);
		for (k = 0; k < p->nei[j].n; ++k) {
			kputl(p->nei[j].a[k].x, out); kputc(',', out); kputw((int32_t)p->nei[j].a[k].y, out);
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

mog_t *mog_read(const char *fn, const mogopt_t *opt)
{
	gzFile fp;
	kseq_t *seq;
	int64_t tot_len = 0, n_arcs = 0, n_tips = 0, n_arc_drop = 0;
	mog128_v nei;
	mog_t *g;

	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	kv_init(nei);
	g = calloc(1, sizeof(mog_t));
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i, j;
		char *q;
		mognode_t *p;
		kv_pushp(mognode_t, g->nodes, &p);
		kv_init(p->nei[0]); kv_init(p->nei[1]); kv_init(p->mapping);
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
			while (isdigit(*q)) { // parse the neighbors
				mog128_t *r;
				kv_push(mog128_t, nei, &r);
				r->x = strtol(q, &q, 10); ++q;
				r->y = strtol(q, &q, 10); ++q;
				g->min_ovlp = g->min_ovlp < r->y? g->min_ovlp : r->y;
				if (max < r->y) max = max2, max = r->y;
				else if (max2 < r->y) max2 = r->y;
			}
			++q; // skip the tailing blank
			thres = (int)(max2 * opt->diff_ratio + .499);
			for (i = 0; i < nei.n; ++i)
				if (nei.a[i].y < thres) nei.a[i].y = 0; // to be deleted in rmdup_128v()
			v128_rmdup(&nei);
			v128_cap(&nei, opt->max_arc);
			kv_copy(mog128_t, p->nei[j], nei);
		}
		// test if to cut a tip
		p->len = seq->seq.l;
		if (opt->flag & MOG_F_DROP_TIP) {
			if ((p->nei[0].n & p->nei[1].n) && p->len < opt->min_el && p->nsr == 1) {
				free(p->nei[0].a); free(p->nei[1].a); // only ->nei[2] have been allocated so far
				--g->nodes.n;
				continue;
			}
		}
		// set ->{seq,cov,max_len}
		p->max_len = p->len + 1;
		kroundup32(p->max_len);
		p->seq = malloc(p->max_len);
		p->cov = malloc(p->max_len);
		strcpy(p->seq, seq->seq.s);
		strcpy(p->cov, seq->qual.s);
		p->aux[0] = p->aux[1] = -1;
	}
	// free and finalize the graph
	kseq_destroy(seq);
	gzclose(fp);
	free(nei.a);
	g->h = build_hash(&g->nodes);
	msg_amend(g);
	//g->rdist = fmg_compute_rdist(&g->nodes);
	return g;
}
