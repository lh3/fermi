#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "priv.h"
#include "utils.h"
#include "kstring.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_INIT2(64,, khint64_t, uint64_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

#define fm128_xlt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y > (b).y))
#define fm128_ylt(a, b) ((int64_t)(a).y > (int64_t)(b).y)
#include "ksort.h"
KSORT_INIT(128x, fm128_t, fm128_xlt)
KSORT_INIT(128y, fm128_t, fm128_ylt)

#define MAX_DEBUBBLE_DIFF 3

typedef khash_t(64) hash64_t;

static hash64_t *build_hash(const fmnode_v *nodes);
static void flip(fmnode_t *p, hash64_t *h);
double fmg_compute_rdist(const fmnode_v *n);

#define arc_mark_del(_x) ((_x).x = (uint64_t)-2, (_x).y = 0)
#define arc_is_del(_x)   ((_x).x == (uint64_t)-2 || (_x).y == 0)

void msg_write_node(const fmnode_t *p, long id, kstring_t *out)
{
	int j, k;
	if (p->l <= 0) return;
	kputc('@', out); kputl(id, out); kputc('_', out); kputw(p->n, out); kputc('_', out); kputl(p->l, out);
	for (j = 0; j < 2; ++j) {
		kputc('\t', out);
		kputl(p->k[j], out); kputc('>', out);
		for (k = 0; k < p->nei[j].n; ++k) {
			if (k) kputc(',', out);
			kputl(p->nei[j].a[k].x, out); kputc(':', out); kputw((int32_t)p->nei[j].a[k].y, out);
		}
		if (p->nei[j].n == 0) kputc('.', out);
	}
	if (p->mapping.n) {
		kputc('\t', out);
		for (j = 0; j < p->mapping.n; ++j) {
			fm128_t *q = &p->mapping.a[j];
			if (j) kputc(',', out);
			kputl(q->x, out); kputc(':', out);
			kputc((q->y&1)? '-' : '+', out); kputc(':', out);
			kputw(q->y>>32, out); kputc(':', out);
			kputw(q->y<<32>>33, out);
		}
	}
	kputc('\n', out);
	ks_resize(out, out->l + 2 * p->l + 5);
	for (j = 0; j < p->l; ++j) out->s[out->l++] = "ACGT"[(int)p->seq[j] - 1];
	out->s[out->l] = 0;
	kputsn("\n+\n", 3, out);
	kputsn(p->cov, p->l, out);
	kputc('\n', out);
}

void msg_print(const fmnode_v *nodes)
{
	size_t i;
	kstring_t out;
	out.l = out.m = 0; out.s = 0;
	for (i = 0; i < nodes->n; ++i) {
		if (nodes->a[i].l > 0) {
			out.l = 0;
			if (out.m) out.s[0] = 0;
			msg_write_node(&nodes->a[i], i, &out);
			fputs(out.s, stdout);
		}
	}
	free(out.s);
}

void msg_nodecpy(fmnode_t *dst, const fmnode_t *src)
{
	int i, max_l;
	memcpy(dst, src, sizeof(fmnode_t));
	max_l = src->l + 1;
	kroundup32(max_l); // this is required; otherwise merge will fail
	dst->seq = calloc(max_l, 1);
	dst->cov = calloc(max_l, 1);
	memcpy(dst->seq, src->seq, src->l);
	memcpy(dst->cov, src->cov, src->l);
	kv_init(dst->nei[0]); kv_init(dst->nei[1]); kv_init(dst->mapping);
	for (i = 0; i < src->nei[0].n; ++i) kv_push(fm128_t, dst->nei[0], src->nei[0].a[i]);
	for (i = 0; i < src->nei[1].n; ++i) kv_push(fm128_t, dst->nei[1], src->nei[1].a[i]);
	for (i = 0; i < src->mapping.n; ++i) kv_push(fm128_t, dst->mapping, src->mapping.a[i]);
}

static inline void rmdup_128v(fm128_v *r)
{
	int l, cnt;
	uint64_t x;
	if (r->n < 2) return;
	ks_introsort(128x, r->n, r->a);
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
	if (cnt) { // we need to shrink the array
		for (l = 0, cnt = 0; l < r->n; ++l)
			if (!arc_is_del(r->a[l])) r->a[cnt++] = r->a[l];
		r->n = cnt;
	}
}

static void rm_dup_arc(fmnode_t *p)
{
	rmdup_128v(&p->nei[0]);
	rmdup_128v(&p->nei[1]);
}

static hash64_t *build_hash(const fmnode_v *nodes)
{
	size_t i;
	int j, ret;
	hash64_t *h;
	double tcpu = cputime();
	h = kh_init(64);
	for (i = 0; i < nodes->n; ++i) {
		const fmnode_t *p = &nodes->a[i];
		if (p->l < 0) continue;
		for (j = 0; j < 2; ++j) {
			khint_t k = kh_put(64, h, p->k[j], &ret);
			if (ret == 0) {
				if (fm_verbose >= 2)
					fprintf(stderr, "[W::%s] tip %ld is duplicated.\n", __func__, (long)p->k[j]);
				kh_val(h, k) = (uint64_t)-1;
			} else kh_val(h, k) = i<<1|j;
		}
	}
	if (fm_verbose >= 3)
		fprintf(stderr, "[%s] build the hash table in %.2f sec\n", __func__, cputime() - tcpu);
	return h;
}

msg_t *msg_read(const char *fn, int drop_tip, int max_arc, float diff_ratio)
{
	gzFile fp;
	kseq_t *seq;
	int64_t tot_len = 0, n_arcs = 0, n_tips = 0, n_arc_drop = 0;
	double tcpu = cputime();
	fm128_v nei;
	msg_t *g;

	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	kv_init(nei);
	g = calloc(1, sizeof(msg_t));
	g->min_ovlp = 1<<31;
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		fmnode_t *p;
		int i, j;
		uint64_t sum;
		uint32_t tmp;
		char *q;
		kv_pushp(fmnode_t, g->nodes, &p);
		kv_init(p->nei[0]); kv_init(p->nei[1]); kv_init(p->mapping);
		if ((q = strchr(seq->name.s, '_')) != NULL) p->n = atoi(q + 1);
		p->l = seq->seq.l;
		tmp = p->l + 1;
		kroundup32(tmp);
		p->seq = malloc(tmp);
		for (i = 0; i < p->l; ++i) p->seq[i] = seq_nt6_table[(int)seq->seq.s[i]];
		p->cov = malloc(tmp);
		strcpy(p->cov, seq->qual.s);
		for (i = 0, sum = 0; i < p->l; ++i) sum += p->cov[i] - 33;
		p->avg_cov = (double)sum / p->l;
		p->aux[0] = p->aux[1] = -1;
		for (j = 0, q = seq->comment.s; j < 2; ++j) {
			int ori_n, max_ovlp = 0, max_ovlp2 = 0, ovlp_thres;
			p->k[j] = strtol(q, &q, 10);
			nei.n = 0;
			if (*q == '>' && q[1] != '.') {
				fm128_t x;
				do {
					++q;
					x.x = strtol(q, &q, 10); ++q;
					x.y = strtol(q, &q, 10);
					g->min_ovlp = g->min_ovlp < x.y? g->min_ovlp : x.y;
					if (max_ovlp < x.y) max_ovlp2 = max_ovlp, max_ovlp = x.y;
					else if (max_ovlp2 < x.y) max_ovlp2 = x.y;
					kv_push(fm128_t, nei, x);
				} while (*q == ',');
				++q;
			} else q += 2; // no arcs
			n_arcs += (ori_n = nei.n);
			ovlp_thres = (int)(max_ovlp2 * diff_ratio + .499);
			for (i = 0; i < nei.n; ++i)
				if (nei.a[i].y < ovlp_thres) nei.a[i].y = 0; // to be deleted in rmdup_128v()
			rmdup_128v(&nei);
			if (nei.n > max_arc) { // excessive connections; keep the top ones
				int thres;
				ks_introsort(128y, nei.n, nei.a);
				thres = nei.a[max_arc].y;
				for (i = 0; i < nei.n; ++i)
					if (nei.a[i].y == thres) break;
				nei.n = i;
			}
			n_arc_drop += ori_n - nei.n;
			kv_copy(fm128_t, p->nei[j], nei);
		}
		if (seq->comment.l - (q - seq->comment.s) >= 5) { // read mapping
			for (; *q && !isdigit(*q); ++q); // skip non-digits
			p->mapping.n = 0;
			while (isdigit(*q)) {
				fm128_t x;
				x.x = strtol(q, &q, 10); ++q;
				x.y = *q == '-'? 1 : 0; q += 2;
				x.y |= (uint64_t)strtol(q, &q, 10)<<32; ++q;
				x.y |= strtol(q, &q, 10)<<1;
				kv_push(fm128_t, p->mapping, x);
				if (*q++ == 0) break;
			}
		}
		if ((p->nei[0].n == 0 || p->nei[1].n == 0) && p->avg_cov < 1.001) {
			if (drop_tip) {
				free(p->nei[0].a); free(p->nei[1].a); free(p->cov); free(p->seq);
				--g->nodes.n;
			}
			++n_tips;
		} else {
			// finalize
			tot_len += seq->seq.l;
			if (fm_verbose >= 4 && g->nodes.n % 100000 == 0)
				fprintf(stderr, "[%s] read %ld nodes in %ld bp in %.2f sec (#tips: %ld; #arcs: %ld; #dropped: %ld)\n",
						__func__, (long)g->nodes.n, (long)tot_len, cputime() - tcpu, (long)n_tips, (long)n_arcs, (long)n_arc_drop);
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	free(nei.a);
	if (fm_verbose >= 3)
		fprintf(stderr, "[%s] finished reading %ld nodes in %ld bp in %.2f sec (#tips: %ld; #arcs: %ld; #dropped: %ld)\n",
				__func__, (long)g->nodes.n, (long)tot_len, cputime() - tcpu, (long)n_tips, (long)n_arcs, (long)n_arc_drop);
	g->h = build_hash(&g->nodes);
	g->rdist = fmg_compute_rdist(&g->nodes);
	msg_amend(g);
	return g;
}

void msg_destroy(msg_t *g)
{
	size_t i;
	for (i = 0; i < g->nodes.n; ++i) {
		fmnode_t *p = &g->nodes.a[i];
		if (p->l < 0) continue;
		free(p->nei[0].a); free(p->nei[1].a); free(p->mapping.a);
		free(p->seq); free(p->cov);
	}
	free(g->nodes.a);
	kh_destroy(64, g->h);
}

static inline uint64_t get_node_id(hash64_t *h, uint64_t tid)
{
	khint_t k;
	k = kh_get(64, h, tid);
	if (k == kh_end(h)) return (uint64_t)(-1);
	return kh_val(h, k);
}

void msg_amend(msg_t *g)
{
	size_t i;
	int j, l, ll;
	double tcpu;
	tcpu = cputime(); // excluding time spent on building the hash table
	for (i = 0; i < g->nodes.n; ++i) {
		fmnode_t *p = &g->nodes.a[i];
		fm128_v *r;
		if (p->l <= 0) continue;
		for (j = 0; j < 2; ++j) {
			int cnt0 = 0, cnt1 = 0;
			for (l = 0; l < p->nei[j].n; ++l) {
				uint64_t x = p->nei[j].a[l].x;
				uint64_t z = get_node_id(g->h, x);
				if (z == (uint64_t)-1) {
					if (fm_verbose >= 5) fprintf(stderr, "[W::%s] tip %ld is non-existing.\n", __func__, (long)x);
					arc_mark_del(p->nei[j].a[l]); // FIXME: should we set -1 or delete it?
					++cnt0;
					continue;
				}
				r = &g->nodes.a[z>>1].nei[z&1];
				for (ll = 0; ll < r->n; ++ll)
					if (r->a[ll].x == p->k[j]) break;
				if (ll == r->n) {
					if (fm_verbose >= 5) fprintf(stderr, "[W::%s] have %ld->%ld but no reverse.\n", __func__, (long)p->k[j], (long)x);
					p->nei[j].a[l].x = (uint64_t)-1;
					++cnt1;
					continue;
				}
			}
			rmdup_128v(&p->nei[j]);
		}
		if (fm_verbose >= 4 && (i+1) % 100000 == 0)
			fprintf(stderr, "[M::%s] amended %ld nodes in %.2f sec\n", __func__, i+1, cputime() - tcpu);
	}
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] amended the graph in %.2f sec\n", __func__, cputime() - tcpu);
}

/*********************
 * remove nodes/arcs *
 *********************/

static void cut_arc(msg_t *g, uint64_t u, uint64_t v, int remove) // delete v from u
{
	int i, j;
	uint64_t x;
	fm128_v *r;
	fmnode_t *p;
	x = get_node_id(g->h, u);
	if (x == (uint64_t)-1) return;
	p = &g->nodes.a[x>>1];
	r = &p->nei[x&1];
	if (remove) {
		for (j = i = 0; j < r->n; ++j)
			if (r->a[j].x != v) r->a[i++] = r->a[j];
		r->n = i;
	} else {
		for (j = 0; j < r->n; ++j)
			if (r->a[j].x == v) arc_mark_del(r->a[j]);
	}
}

static void drop_arc(msg_t *g, size_t id, int min_ovlp, float min_ovlp_ratio)
{
	fmnode_t *p = &g->nodes.a[id];
	int j, l, cnt;
	if (p->l < 0 || (min_ovlp == 0 && (min_ovlp_ratio <= 0.01 || min_ovlp_ratio >= 0.99))) return;
	for (j = 0; j < 2; ++j) {
		fm128_v *r = &p->nei[j];
		int max = 0;
		if (r->n == 0) continue;
		for (l = 0; l < r->n; ++l)
			if (r->a[l].y > max) max = r->a[l].y;
		for (l = cnt = 0; l < r->n; ++l) {
			if (r->a[l].y < min_ovlp || (double)r->a[l].y/max < min_ovlp_ratio) {
				if (r->a[l].x != p->k[0] && r->a[l].x != p->k[1])
					cut_arc(g, r->a[l].x, p->k[j], 1);
				arc_mark_del(r->a[l]);
				++cnt;
			}
		}
		if (cnt) {
			for (l = cnt = 0; l < r->n; ++l)
				if (!arc_is_del(r->a[l])) r->a[cnt++] = r->a[l];
			r->n = cnt;
		}
	}
}

static void drop_all_weak_arcs(msg_t *g, int min_ovlp, float min_ovlp_ratio)
{
	size_t i;
	double t = cputime();
	for (i = 0; i < g->nodes.n; ++i)
		drop_arc(g, i, min_ovlp, min_ovlp_ratio);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] cut weak arcs with threshold (%d,%.2f) in %.3f sec.\n",
				__func__, min_ovlp, min_ovlp_ratio, cputime() - t);
}

static inline void rmnode_force(msg_t *g, fmnode_t *p)
{
	int i;
	if (p->l < 0) return;
	for (i = 0; i < p->nei[0].n; ++i)
		if (p->nei[0].a[i].x != p->k[0] && p->nei[0].a[i].x != p->k[1])
			cut_arc(g, p->nei[0].a[i].x, p->k[0], 1);
	for (i = 0; i < p->nei[1].n; ++i)
		if (p->nei[1].a[i].x != p->k[0] && p->nei[1].a[i].x != p->k[1])
			cut_arc(g, p->nei[1].a[i].x, p->k[1], 1);
	p->nei[0].n = p->nei[1].n = 0;
	p->l = -1; p->n = 0;
}

static inline void add_arc(msg_t *g, uint64_t u, uint64_t v, int ovlp)
{
	fm128_v *r;
	fm128_t z;
	uint64_t x = get_node_id(g->h, u);
	if (x == (uint64_t)-1) return;
	r = &g->nodes.a[x>>1].nei[x&1];
	z.x = v; z.y = ovlp;
	kv_push(fm128_t, *r, z);
}

static void rmnode_int(msg_t *g, size_t id)
{
	int i, j;
	fmnode_t *p = &g->nodes.a[id];
	if (p->l < 0) return;
	if (p->nei[0].n && p->nei[1].n) {
		for (i = 0; i < p->nei[0].n; ++i) {
			if (p->nei[0].a[i].x == p->k[0] || p->nei[0].a[i].x == p->k[1]) continue;
			for (j = 0; j < p->nei[1].n; ++j) {
				int ovlp = (int)(p->nei[0].a[i].y + p->nei[1].a[j].y) - p->l;
				if (p->nei[1].a[j].x == p->k[0] || p->nei[1].a[j].x == p->k[1]) continue;
				if (ovlp >= g->min_ovlp) {
					add_arc(g, p->nei[0].a[i].x, p->nei[1].a[j].x, ovlp);
					add_arc(g, p->nei[1].a[j].x, p->nei[0].a[i].x, ovlp);
				}
			}
		}
	}
	rmnode_force(g, p);
}

void msg_rm_tips(msg_t *g, int min_len, int min_cnt)
{
	size_t i;
	double t = cputime();
	if (min_len < 0) return;
	if (g->h == 0) g->h = build_hash(&g->nodes);
	for (i = 0; i < g->nodes.n; ++i) {
		fmnode_t *p = &g->nodes.a[i];
		if (p->nei[0].n && p->nei[1].n) continue; // not a tip
		if (p->l < min_len && p->n < min_cnt) rmnode_force(g, p);
	}
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] removed tips shorter than %dbp in %.3f sec.\n", __func__, min_len, cputime() - t);
}

/******************
 * Bubble popping *
 ******************/

typedef struct {
	fm64_v visited, order;
	fm128_v heap;
	hash64_t *h, *kept;
} popaux_t;
/*
static void stretch_simple_circle(fmnode_v *nodes, hash64_t *h, int max_ll, size_t id)
{
	fmnode_t *q, *p = &nodes->a[id];
	fm128_v *r;
	int j, cnt;
	if (p->l < 0) return;
	r = &p->nei[0];
	for (j = cnt = 0; j < r->n; ++j)
		if (r->a[j].x == p->k[1])
			++cnt, r->a[j].x = 0;
	if (cnt) {
		for (j = cnt = 0; j < r->n; ++j)
			if (r->a[j].x) r->a[cnt++] = r->a[j];
		r->n = cnt;
		r = &p->nei[1];
		for (j = cnt = 0; j < r->n; ++j)
			if (r->a[j].x != p->k[0]) r->a[cnt++] = r->a[j];
		assert(cnt + 1 == r->n);
		r->n = cnt;
	}
	if (p->nei[0].n == 1 && p->nei[1].n == 1) { // test a circle like a->b->a (test b)
		uint64_t x[2];
		x[0] = get_node_id(h, p->nei[0].a[0].x);
		x[1] = get_node_id(h, p->nei[1].a[0].x);
		if (x[0] == (uint64_t)-1 || x[1] == (uint64_t)-1) return;
		if ((x[0] ^ x[1]) == 1) {
			int i, new_l, max_l, o[2];
			char *seq, *cov;
			q = &nodes->a[x[0]>>1];
			if (q->avg_cov < 4. || q->l > max_ll || q->l < 0) return; 
			o[x[0]&1] = p->nei[0].a[0].y;
			o[x[1]&1] = p->nei[1].a[0].y;
			if (x[1]&1) flip(p, h);
			new_l = p->l + q->l * 2 - o[0] - o[1];
			assert(new_l > 0);
			max_l = new_l + 1;
			kroundup32(max_l);
			seq = calloc(max_l, 1);
			cov = calloc(max_l, 1);
			memcpy(seq, q->seq, q->l);
			memcpy(seq + (q->l - o[1]), p->seq, p->l);
			memcpy(seq + (q->l + p->l - o[1] - o[0]), q->seq, q->l);
			for (i = 0; i < q->l; ++i) cov[i] = q->cov[i];
			for (; i < new_l; ++i) cov[i] = 33;
			for (i = 0; i < p->l; ++i) {
				int c = cov[i + (q->l - o[1])] - 33 + p->cov[i];
				cov[i + (q->l - o[1])] = c > 126? 126 : c;
			}
			for (i = 0; i < q->l; ++i) {
				int c = cov[i + p->l + q->l - o[0] - o[1]] - 33 + q->cov[i];
				cov[i + p->l + q->l - o[0] - o[1]] = c > 126? 126 : c;
			}
			rmnode(nodes, h, id);
			free(q->seq); free(q->cov);
			//fprintf(stderr, "%d, %d, %d, %d\n", o[0], o[1], q->l, new_l);
			q->l = new_l;
			q->seq = seq; q->cov = cov;
			q->seq[new_l] = q->cov[new_l] = 0;
		} else if (p->l < max_ll) {
			cut_arc(nodes, h, p->nei[0].a[0].x, p->nei[1].a[0].x, 1);
			cut_arc(nodes, h, p->nei[0].a[1].x, p->nei[0].a[0].x, 1);
		}
	}
}
*/
static void pop_complex_bubble(msg_t *g, size_t idd, int max_len, int max_nodes, popaux_t *aux)
{
	fmnode_t *q, *p = &g->nodes.a[idd>>1];
	fm128_t u, v;
	fm128_v *r;
	khint_t k;
	uint64_t tmp, term = (uint64_t)-1;
	int i, j, ret, has_loop;

	if (p->l < 0 || p->nei[idd&1].n < 2) return;
	// initialize the aux structure
	aux->heap.n = aux->visited.n = aux->order.n = 0;
	kh_clear(64, aux->h);
	p->aux[idd&1] = 0;
	v.x = idd; v.y = 0;
	kv_push(fm128_t, aux->heap, v);
	kh_put(64, aux->h, p->k[idd&1], &ret);
	// Dijkstra-like graph traversal
	while (aux->heap.n) {
		if (aux->heap.n > max_nodes) goto end_popcomp;
		// take the shortest extension from the top of the heap
		u = aux->heap.a[0];
		p = &g->nodes.a[u.x>>1];
		//fprintf(stderr, "%lld,%lld; [%lld,%lld]\n", u.x, u.y, nodes->a[u.x>>1].k[0], nodes->a[u.x>>1].k[1]);
		aux->heap.a[0] = aux->heap.a[--aux->heap.n];
		ks_heapdown_128y(0, aux->heap.n, aux->heap.a);

		r = &p->nei[u.x&1];
		if (r->n > max_nodes) goto end_popcomp;
		if (term == (uint64_t)-1 && aux->heap.n == 0 && u.x != idd) {
			term = u.x;
			break;
		} else if (r->n == 0 || p->aux[u.x&1] > max_len) { // hit a dead-end or exceed the length limit
			if (term != (uint64_t)-1) {
				if (term != u.x) goto end_popcomp; // multiple dead-end; stop
			} else term = u.x; // the first dead-end
		} else { // try to extend
			kv_push(uint64_t, aux->order, u.x);
			for (i = 0; i < r->n; ++i) {
				v.x = get_node_id(g->h, r->a[i].x);
				if (v.x == (uint64_t)-1) continue;
				v.x ^= 1;
				q = &g->nodes.a[v.x>>1];
				//fprintf(stderr, "[%c] %lld[%lld,%lld] -> %lld[%lld,%lld]\n", "-+"[q->aux[v.x&1]<0], u.x,
				//	nodes->a[u.x>>1].k[0], nodes->a[u.x>>1].k[1], v.x^1, nodes->a[v.x>>1].k[0], nodes->a[v.x>>1].k[1]);
				if (q->aux[v.x&1] < 0) { // has not been visited
					v.y = q->aux[v.x&1] = p->aux[u.x&1] - r->a[i].y + q->l;
					kv_push(fm128_t, aux->heap, v);
					ks_heapup_128y(aux->heap.n, aux->heap.a);
					kv_push(uint64_t, aux->visited, v.x);
					kh_put(64, aux->h, q->k[0], &ret);
					kh_put(64, aux->h, q->k[1], &ret);
				}
			}
		}
	}
	// check if this is really a bubble
	for (i = 0; i < aux->visited.n; ++i) {
		if (aux->visited.a[i] == idd) continue;
		p = &g->nodes.a[aux->visited.a[i]>>1];
		r = &p->nei[(aux->visited.a[i]^1)&1];
		for (j = 0; j < r->n; ++j) {
			khint_t k = kh_get(64, aux->h, r->a[j].x);
			if (k == kh_end(aux->h)) goto end_popcomp;
		}
	}
	// prepare for finding the longest path
	for (i = 0; i < aux->visited.n; ++i) {
		p = &g->nodes.a[aux->visited.a[i]>>1];
		p->aux[0] = p->aux[1] = 0;
	}
	kh_clear(64, aux->h);
	kh_clear(64, aux->kept);
	// find the longest path. FIXME: this is an approximate implementation and may FAIL!
	for (i = 0, has_loop = 0; i < aux->order.n; ++i) {
		uint64_t x = aux->order.a[i];
		kh_put(64, aux->kept, x>>1, &ret); // here aux->kept is used to break loop(s)
		p = &g->nodes.a[x>>1];
		r = &p->nei[x&1];
		for (j = 0; j < r->n; ++j) {
			uint64_t y = get_node_id(g->h, r->a[j].x);
			int d;
			if (y == (uint64_t)-1) continue;
			y ^= 1;
			q = &g->nodes.a[y>>1];
			d = (int)(q->avg_cov * q->l + .499);
			if (q->aux[y&1] < p->aux[x&1] + d) {
				q->aux[y&1] = p->aux[x&1] + d;
				if (kh_get(64, aux->kept, y>>1) == kh_end(aux->kept)) { // do not set if there is a loop
					k = kh_put(64, aux->h, y^1, &ret);
					kh_val(aux->h, k) = x;
				} else has_loop = 1;
			}
		}
	}
	if (has_loop) goto end_popcomp;
	// backtrack
	aux->heap.n = 0;
	tmp = term;
	kh_clear(64, aux->kept);
	while (1) {
		kh_put(64, aux->kept, tmp>>1, &ret);
		k = kh_get(64, aux->h, tmp^1);
		if (k == kh_end(aux->h)) break;
		tmp = kh_val(aux->h, k);
	}
	// delete nodes
	for (j = 0; j < aux->visited.n; ++j) {
		khint_t k = kh_get(64, aux->kept, aux->visited.a[j]>>1);
		if (k == kh_end(aux->kept)) {
			rmnode_force(g, &g->nodes.a[aux->visited.a[j]>>1]);
			//fprintf(stderr, "del %lld[%lld,%lld]\n", aux->visited.a[j], nodes->a[aux->visited.a[j]>>1].k[0], nodes->a[aux->visited.a[j]>>1].k[1]);
		}
	}
end_popcomp:
	for (i = 0; i < aux->visited.n; ++i) {
		p = &g->nodes.a[aux->visited.a[i]>>1];
		p->aux[0] = p->aux[1] = -1;
	}
}

void msg_popbub(msg_t *g, int max_len, int max_nodes)
{
	fm128_v tmp;
	size_t i;
	double t = cputime();
	popaux_t *aux;

	kv_init(tmp);
	aux = calloc(1, sizeof(popaux_t));
	aux->h = kh_init(64);
	aux->kept = kh_init(64);
	for (i = 0; i < g->nodes.n; ++i) {
		fm128_t x;
		fmnode_t *p;
		//stretch_simple_circle(nodes, h, max_len, i);
		p = &g->nodes.a[i];
		if (p->l > 0 && (p->nei[0].n > 1 || p->nei[1].n > 1)) {
			x.x = i; x.y = g->nodes.a[i].l;
			kv_push(fm128_t, tmp, x);
		}
	}
	ks_introsort(128y, tmp.n, tmp.a);
	for (i = 0; i < tmp.n; ++i) {
		pop_complex_bubble(g, tmp.a[i].x<<1|0, max_len, max_nodes, aux);
		pop_complex_bubble(g, tmp.a[i].x<<1|1, max_len, max_nodes, aux);
	}

	free(tmp.a);
	kh_destroy(64, aux->h);
	kh_destroy(64, aux->kept);
	free(aux->heap.a); free(aux->visited.a);
	if (fm_verbose >= 3)
		fprintf(stderr, "[M::%s] popped complex bubbles in %.3f sec.\n", __func__, cputime() - t);
}

static void pop_simple_bubble(msg_t *g, size_t id, double min_bub_ratio, double min_bub_cov)
{
	fmnode_t *p = &g->nodes.a[id], *q[2], *top_p, *tmp_p;
	fm128_v *r[2];
	uint64_t nei[2];
	int j, cnt;
	double max_cov;
	uint64_t top_id = (uint64_t)-1;

	if (p->l <= 0 || p->nei[0].n != 1 || p->nei[1].n != 1) return;
	for (j = 0; j < 2; ++j) {
		nei[j] = get_node_id(g->h, p->nei[j].a[0].x);
		if (nei[j] == (uint64_t)-1) return;
		q[j] = &g->nodes.a[nei[j]>>1];
	}
	if (q[0]->l < 0 || q[1]->l < 0) return; // deleted node
	if (q[0]->nei[nei[0]&1].n <= 1 || q[1]->nei[nei[1]&1].n <= 1) return; // unmerged or inconsistent arc
	r[0] = &q[0]->nei[nei[0]&1]; r[1] = &q[1]->nei[nei[1]&1];
	for (j = 0, max_cov = 0, cnt = 0; j < r[0]->n; ++j) {
		uint64_t t = get_node_id(g->h, r[0]->a[j].x);
		if (t == (uint64_t)-1) continue;
		tmp_p = &g->nodes.a[t>>1];
		if (tmp_p->nei[0].n != 1 || tmp_p->nei[1].n != 1) continue; // skip this node
		if (t&1) flip(tmp_p, g->h); // s.t. tmp_p is on the same strand as p
		if (tmp_p->nei[1].a[0].x != p->nei[1].a[0].x) continue; // not a multi-edge
		if (tmp_p->avg_cov > max_cov) max_cov = tmp_p->avg_cov, top_id = t;
		++cnt;
	}
	if (cnt < 2) return;
	assert(top_id != (uint64_t)-1);
	top_p = &g->nodes.a[top_id>>1];
	for (j = 0; j < r[0]->n; ++j) {
		uint64_t t = get_node_id(g->h, r[0]->a[j].x);
		int l, diff, ml = 0, to_del = 0, beg[2], end[2];
		double cov[2];
		if (t == (uint64_t)-1) continue;
		if (top_id == t) continue; // we do not process the node with the highest coverage
		tmp_p = &g->nodes.a[t>>1];
		if (tmp_p->nei[0].n != 1 || tmp_p->nei[1].n != 1) continue; // skip this node
		if (tmp_p->nei[1].a[0].x != p->nei[1].a[0].x) continue; // not a multi-edge
		// the following is really nasty. A banded global alignment would look much cleaner and work better.
		beg[0] = top_p->nei[0].a[0].y; end[0] = top_p->l - top_p->nei[1].a[0].y;
		beg[1] = tmp_p->nei[0].a[0].y; end[1] = tmp_p->l - tmp_p->nei[1].a[0].y;
		for (l = diff = 0; l < end[0] - beg[0] && l < end[1] - beg[1]; ++l)
			if (top_p->seq[l + beg[0]] != tmp_p->seq[l + beg[1]])
				if (diff++ == 0) ml = l;
		if (diff) { // then compute the number of matching bases from the end of the sequence
			for (l = 0; l < end[0] - beg[0] && l < end[1] - beg[1]; ++l)
				if (top_p->seq[end[0] - 1 - l] != tmp_p->seq[end[1] - 1 - l]) break;
			ml += l;
		} else ml = l;
		for (l = beg[0], cov[0] = 0; l < end[0]; ++l) cov[0] += top_p->cov[l] - 33;
		cov[0] /= end[0] > beg[0]? end[0] - beg[0] : 1;
		for (l = beg[1], cov[1] = 0; l < end[1]; ++l) cov[1] += tmp_p->cov[l] - 33;
		cov[1] /= end[1] > beg[1]? end[1] - beg[1] : 1;
		if (top_p != tmp_p && abs((end[0]-beg[0]) - (end[1]-beg[1])) < MAX_DEBUBBLE_DIFF) {
			if (diff < MAX_DEBUBBLE_DIFF) to_del = 1;
			else if (end[0]-beg[0] - ml < MAX_DEBUBBLE_DIFF || end[1]-beg[1] - ml < MAX_DEBUBBLE_DIFF) to_del = 1;
		}
		if (cov[0] > 0. && cov[1] > 0. && (cov[1] / cov[0] >= min_bub_ratio || cov[1] >= min_bub_cov)) to_del = 0;
		if (1&&to_del) {
			cut_arc(g, tmp_p->k[0], q[0]->k[nei[0]&1], 1);
			cut_arc(g, q[0]->k[nei[0]&1], tmp_p->k[0], 0);
			cut_arc(g, tmp_p->k[1], q[0]->k[nei[1]&1], 1);
			cut_arc(g, q[1]->k[nei[1]&1], tmp_p->k[1], 1);
			tmp_p->l = -1;
		}
	}
	for (j = 0, cnt = 0; j < r[0]->n; ++j)
		if (!arc_is_del(r[0]->a[j])) r[0]->a[cnt++] = r[0]->a[j];
	r[0]->n = cnt;
}

/*********************
 * Unambiguous merge *
 *********************/

#define __swap(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))

static void flip(fmnode_t *p, hash64_t *h)
{
	fm128_v t;
	khint_t k;
	seq_revcomp6(p->l, (uint8_t*)p->seq);
	seq_reverse(p->l, (uint8_t*)p->cov);
	__swap(p->k[0], p->k[1]);
	t = p->nei[0]; p->nei[0] = p->nei[1]; p->nei[1] = t;
	k = kh_get(64, h, p->k[0]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
	k = kh_get(64, h, p->k[1]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
	if (p->mapping.n) { // flip the mapping
		int i;
		for (i = 0; i < p->mapping.n; ++i) {
			int beg, end;
			beg = p->mapping.a[i].y>>32;
			end = p->mapping.a[i].y<<32>>33;
			p->mapping.a[i].y = (uint64_t)(p->l - end)<<32 | (p->l - beg)<<1 | ((p->mapping.a[i].y&1)^1);
		}
	}
}

static int merge(msg_t *g, size_t w) // merge i's neighbor to the right-end of i
{
	fmnode_t *p = &g->nodes.a[w], *q;
	khint_t kp, kq;
	uint32_t tmp;
	int i, j, new_l;
	if (p->nei[1].n != 1) return -1; // cannot be merged
	kq = kh_get(64, g->h, p->nei[1].a[0].x);
	if (kq == kh_end((hash64_t*)g->h)) return -2; // not found
	q = &g->nodes.a[kh_val((hash64_t*)g->h, kq)>>1];
	if (p == q) return -4; // this is a loop
	if (q->nei[kh_val((hash64_t*)g->h, kq)&1].n != 1) return -3; // cannot be merged
	// we can perform a merge
	if (kh_val((hash64_t*)g->h, kq)&1) flip(q, g->h); // a "><" bidirectional arc; flip q
	kp = kh_get(64, g->h, p->k[1]); assert(kp != kh_end((hash64_t*)g->h)); // get the iterator to p
	kh_del(64, g->h, kp); kh_del(64, g->h, kq); // remove the two ends of the arc in the hash table
	if (p->k[1] != q->nei[0].a[0].x || q->k[0] != p->nei[1].a[0].x) {
		fprintf(stderr, "[W::%s] a rare bug prevents merging: %ld != %ld or %ld != %ld. Continue anyway.\n",
				__func__, (long)p->k[1], (long)q->nei[0].a[0].x, (long)q->k[0], (long)p->nei[1].a[0].x);
		return -4;
	}
	assert(p->k[1] == q->nei[0].a[0].x);
	assert(q->k[0] == p->nei[1].a[0].x);
	assert(p->nei[1].a[0].y == q->nei[0].a[0].y);
	assert(p->l >= p->nei[1].a[0].y && q->l >= p->nei[1].a[0].y); // "==" may happen due to end trimming
	p->n += q->n;
	new_l = p->l + q->l - p->nei[1].a[0].y;
	tmp = p->l;
	kroundup32(tmp);
	if (new_l + 1 >= tmp) { // then double p->seq and p->cov
		tmp = new_l + 1;
		kroundup32(tmp);
		p->seq = realloc(p->seq, tmp);
		p->cov = realloc(p->cov, tmp);
	}
	// merge seq and cov
	for (i = p->l - p->nei[1].a[0].y, j = 0; j < q->l; ++i, ++j) { // write seq and cov
		p->seq[i] = q->seq[j];
		if (i < p->l) {
			if ((int)p->cov[i] + (q->cov[j] - 33) > 126) p->cov[i] = 126;
			else p->cov[i] += q->cov[j] - 33;
		} else p->cov[i] = q->cov[j];
	}
	p->seq[new_l] = p->cov[new_l] = 0;
	p->avg_cov = (p->l * p->avg_cov + q->l * q->avg_cov) / new_l; // recalculate coverage
	p->l = new_l;
	// merge neighbors
	free(p->nei[1].a);
	p->nei[1] = q->nei[1]; p->k[1] = q->k[1];
	if (p->mapping.n || q->mapping.n) {
		for (i = 0; i < q->mapping.n; ++i)
			kv_push(fm128_t, p->mapping, q->mapping.a[i]);
		ks_introsort_128x(p->mapping.n, p->mapping.a);
		for (i = 1; i < p->mapping.n; ++i)
			if ((p->mapping.a[i].x ^ p->mapping.a[i-1].x) == 1)
				p->mapping.a[i].x = p->mapping.a[i-1].x = 0xffffffffU;
		for (i = j = 0; i < p->mapping.n; ++i)
			if (p->mapping.a[i].x != 0xffffffffU)
				p->mapping.a[j++] = p->mapping.a[i];
		p->mapping.n = j;
	}
	// update the hash table
	kp = kh_get(64, g->h, p->k[1]);
	assert(kp != kh_end((hash64_t*)g->h));
	kh_val((hash64_t*)g->h, kp) = w<<1|1;
	// clean up q
	free(q->cov); free(q->seq);
	q->cov = q->seq = 0; q->l = -1; q->n = 0;
	free(q->nei[0].a); free(q->mapping.a);
	q->nei[0].n = q->nei[0].m = q->nei[1].n = q->nei[1].m = q->mapping.n = q->mapping.m = 0;
	q->nei[0].a = q->nei[1].a = q->mapping.a = 0; // q->nei[1] has been copied over to p->nei[1], so we can delete it
	return 0;
}

void msg_join_unambi(msg_t *g)
{
	size_t i;
	fmnode_v *nodes = &g->nodes;
	double tcpu = cputime();
	for (i = 0; i < nodes->n; ++i) rm_dup_arc(&nodes->a[i]);
	for (i = 0; i < nodes->n; ++i) {
		if (nodes->a[i].l <= 0) continue;
		while (merge(g, i) == 0);
		flip(&g->nodes.a[i], g->h);
		while (merge(g, i) == 0);
	}
	if (fm_verbose >= 2)
		fprintf(stderr, "[M::%s] joined unambiguous arcs in %.2f sec\n", __func__, cputime() - tcpu);
}

/*************************
 * SW high-coverage tips *
 *************************/

#include "ksw.h"

static void tip_sw(msg_t *g, size_t id, int tip_len, int min_cnt)
{
	int i, j, k, dir, max_l, l_qry;
	fmnode_t *p, *q, *t;
	fm128_v *r;
	uint64_t w;
	uint8_t *seq;
	int8_t mat[16];
	ksw_query_t *qry;
	ksw_aux_t aux;

	p = &g->nodes.a[id];
	if (p->l < 0 || p->n < min_cnt || p->l > tip_len) return;
	if (p->nei[0].n + p->nei[1].n != 1) return;
	dir = p->nei[0].n? 0 : 1;
	w = get_node_id(g->h, p->nei[dir].a[0].x);
	if (w == (uint64_t)-1) return;
	q = &g->nodes.a[w>>1];
	if (q->nei[w&1].n == 1) return; // this is already merge-able

	// get the query ready
	max_l = (p->l - p->nei[dir].a[0].y) * 2;
	for (i = k = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? 1 : -3;
	seq = malloc(max_l + 1);
	if (dir == 0) { // forward strand
		for (j = p->nei[dir].a[0].y, k = 0; j < p->l; ++j)
			seq[k++] = p->seq[j] - 1;
	} else { // reverse
		for (j = p->l - p->nei[dir].a[0].y - 1, k = 0; j >= 0; --j)
			seq[k++] = 4 - p->seq[j];
	}
	l_qry = k;
	qry = ksw_qinit(2, l_qry, seq, 4, mat);
	aux.gapo = 3; aux.gape = 2; aux.T = l_qry / 2;
	//fprintf(stderr, "===> %lld:%lld, %d, %ld <===\n", p->k[0], p->k[1], p->n, q->nei[w&1].n);
	//for (j = 0; j < k; ++j) fputc("ACGTN"[(int)seq[j]], stderr); fputc('\n', stderr);

	r = &q->nei[w&1];
	for (i = 0; i < r->n; ++i) {
		if (r->a[i].x == p->k[dir]) continue;
		w = get_node_id(g->h, r->a[i].x);
		if (w == (uint64_t)-1) continue;
		// get the target sequence
		t = &g->nodes.a[w>>1];
		if (w&1) { // reverse strand
			for (j = t->l - r->a[i].y - 1, k = 0; j >= 0 && k < max_l; --j)
				seq[k++] = 4 - t->seq[j];
		} else {
			for (j = r->a[i].y, k = 0; j < t->l && k < max_l; ++j)
				seq[k++] = t->seq[j] - 1;
		}
		ksw_sse2(qry, k, seq, &aux);
		//for (j = 0; j < k; ++j) fputc("ACGTN"[(int)seq[j]], stderr); fprintf(stderr, "\t%d\n", aux.score);
		if (aux.score) {
			double r_diff, n_diff;
			n_diff = .25 * (l_qry - aux.score);
			r_diff = n_diff / l_qry;
			if ((int)(n_diff + .499) == 1 || r_diff < 0.1) break;
		}
	}
	if (i != r->n) rmnode_force(g, p);

	free(seq); free(qry);
}

/*********************************************
 * A-statistics and simplistic flow analysis *
 *********************************************/

#define A_THRES 20.
#define A_MIN_SUPP 5

double fmg_compute_rdist(const fmnode_v *n)
{
	int j;
	uint64_t *srt;
	double rdist = -1., t;
	int64_t i, sum_n_all, sum_n, sum_l;

	t = cputime();
	srt = calloc(n->n, 8);
	for (i = 0, sum_n_all = 0; i < n->n; ++i) {
		srt[i] = (uint64_t)n->a[i].n<<32 | i;
		sum_n_all += n->a[i].n;
	}
	ks_introsort_uint64_t(n->n, srt);

	for (j = 0; j < 2; ++j) {
		sum_n = sum_l = 0;
		for (i = n->n - 1; i >= 0; --i) {
			const fmnode_t *p = &n->a[srt[i]<<32>>32];
			int tmp1, tmp2;
			tmp1 = tmp2 = 0;
			if (p->nei[0].n) ++tmp1, tmp2 += p->nei[0].a[0].y;
			if (p->nei[1].n) ++tmp1, tmp2 += p->nei[1].a[0].y;
			if (tmp1) tmp2 /= tmp1;
			if (rdist > 0.) {
				double A = (p->l - tmp1) / rdist - p->n * M_LN2;
				if (A < A_THRES) continue;
			}
			sum_n += p->n;
			sum_l += p->l - tmp1;
			if (sum_n >= sum_n_all * 0.75) break;
		}
		rdist = (double)sum_l / sum_n;
	}
	if (fm_verbose >= 3) {
		fprintf(stderr, "[M::%s] average read distance %.3f, computed in %.3f seconds.\n", __func__, rdist, cputime() - t);
		fprintf(stderr, "[M::%s] approximate genome size: %.0f (likely an underestimate)\n", __func__, rdist * sum_n_all);
	}

	free(srt);
	return rdist;
}

static void flow_flt1(fmgraph_t *g, size_t idd, double thres)
{
	fmnode_t *p, *q, *t;
	double A;
	uint64_t u, v;
	fm128_v *r;
	int i, n;

	p = &g->nodes.a[idd>>1];
	if (p->nei[idd&1].n != 1) return; // deg != 1
	A = (p->l - p->nei[idd&1].a[0].y) / g->rdist - p->n * M_LN2;
	if (A < thres && p->n < A_MIN_SUPP) return; // p not significantly unique
	u = get_node_id(g->h, p->nei[idd&1].a[0].x);
	if (u == (uint64_t)-1) return; // internal error
	q = &g->nodes.a[u>>1];
	if (p == q) return;
	if (q->nei[u&1].n < 2) return; // well, p and q can be merged already
	assert(q->l >= p->nei[idd&1].a[0].y);
	A = (q->l - p->nei[idd&1].a[0].y) / g->rdist - q->n * M_LN2;
	if (A < thres) return; // q not significantly unique
//	fprintf(stderr, "%lld:%lld, %lld:%lld\n", p->k[0], p->k[1], q->k[0], q->k[1]);

	r = &q->nei[u&1];
	for (i = n = 0; i < r->n; ++i) {
		int to_cut = 0;
		v = get_node_id(g->h, r->a[i].x);
		if (v == (uint64_t)-1) continue;
		t = &g->nodes.a[v>>1];
		if (t->nei[v&1].n == 1) {
			// Should we also consider to cut in this case? Perhaps does not matter.
			if (t->n <= 2 && t->nei[(v&1)^1].n == 0) to_cut = 1; // a small tip
		} else to_cut = 1;
		if (to_cut) {
			if (r->a[i].x != q->k[0] && r->a[i].x != q->k[1])
				cut_arc(g, r->a[i].x, q->k[u&1], 1); // remove q from its neighbor
			arc_mark_del(r->a[i]);
			++n;
		}
	}
	if (n) {
		for (i = n = 0; i < r->n; ++i)
			if (!arc_is_del(r->a[i]))
				r->a[n++] = r->a[i];
		r->n = n;
	}
}

static void flow_flt(fmgraph_t *g, double thres)
{
	size_t i;
	for (i = 0; i < g->nodes.n; ++i) {
		flow_flt1(g, i<<1|0, thres);
		flow_flt1(g, i<<1|1, thres);
	}
}

/**************
 * Key portal *
 **************/

void msg_clean(msg_t *g, const fmclnopt_t *opt)
{
	fm64_v stack;
	size_t i;
	int j;

	kv_init(stack);
	if (g->h == 0) g->h = build_hash(&g->nodes);
//	if (opt->uniq_cut) { flow_flt(g, A_THRES * 2.); msg_join_unambi(g); }
	for (j = 0; j < opt->n_iter; ++j) {
		double r = opt->n_iter == 1? 1. : .5 + .5 * j / (opt->n_iter - 1);
		drop_all_weak_arcs(g, opt->min_ovlp * r, opt->min_ovlp_ratio * r);
		msg_rm_tips(g, opt->min_ext_len * r, opt->min_ext_cnt * r > 2.? opt->min_ext_cnt * r : 2);
		msg_join_unambi(g);
	}
	if (g->min_ovlp < opt->min_ovlp) g->min_ovlp = opt->min_ovlp;
	for (i = 0; i < g->nodes.n; ++i)
		tip_sw(g, i, opt->min_ext_len, opt->min_ext_cnt);
	if (opt->min_int_cnt >= 2) {
		double t = cputime();
		for (i = 0; i < g->nodes.n; ++i)
			if (g->nodes.a[i].n < opt->min_int_cnt)
				rmnode_int(g, i);
		for (i = 0; i < g->nodes.n; ++i) rm_dup_arc(&g->nodes.a[i]);
		fprintf(stderr, "[M::%s] removed weak arcs in %.3f sec.\n", __func__, cputime() - t);
		drop_all_weak_arcs(g, opt->min_ovlp, opt->min_ovlp_ratio);
		msg_rm_tips(g, opt->min_ext_len, opt->min_ext_cnt);
		msg_join_unambi(g);
	}
	flow_flt(g, opt->a_thres > 20.? opt->a_thres : 20.);
	if (opt->aggressive_pop) {
		msg_popbub(g, opt->max_bub_len, 25);
		msg_rm_tips(g, opt->min_ext_len, opt->min_ext_cnt);
		msg_join_unambi(g);
	} else if (opt->min_bub_cov >= 1. && opt->min_bub_ratio < 1.) {
		double t = cputime();
		for (i = 0; i < g->nodes.n; ++i)
			pop_simple_bubble(g, i, opt->min_bub_ratio, opt->min_bub_cov);
		fprintf(stderr, "[M::%s] popped simple bubbles in %.3f sec.\n", __func__, cputime() - t);
		msg_rm_tips(g, opt->min_ext_len, opt->min_ext_cnt);
		msg_join_unambi(g);
	}
	// free
	free(stack.a);
}
