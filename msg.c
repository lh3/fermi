#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fermi.h"
#include "rld.h"
#include "utils.h"
#include "kstring.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khash.h"
KHASH_MAP_INIT_INT64(64, uint64_t)
KHASH_INIT2(arc,, uint64_t, int, 1, kh_int64_hash_func, kh_int64_hash_equal)

#define fm128_xlt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y > (b).y))
#define fm128_ylt(a, b) ((a).y > (b).y)
#include "ksort.h"
KSORT_INIT(128x, fm128_t, fm128_xlt)
KSORT_INIT(128y, fm128_t, fm128_ylt)

#define MAX_DEBUBBLE_DIFF 3
#define MAX_POP_EXTENSION 1000
#define MAX_NEIGHBORS     512

typedef struct {
	uint64_t k[2];
	kh_arc_t *arc[2];
	int l, aux[2];
	float avg_cov;
	char *seq, *cov;
} fmnode_t, *fmnode_p;

typedef khash_t(64) hash64_t;

typedef struct {
	kvec_t(fmnode_p) nodes;
	hash64_t *h;
	int min_ovlp;
} msg_t;

static void flip(fmnode_t *p, hash64_t *h);
static void pop_simple_bubble(msg_t *g, size_t id, double min_bub_ratio, double min_bub_cov, int max_bubble_diff);

fmnode_t *msg_node_init()
{
	fmnode_t *p;
	p = calloc(1, sizeof(fmnode_t));
	p->arc[0] = kh_init(arc);
	p->arc[1] = kh_init(arc);
	return p;
}

void msg_node_destroy(fmnode_t *p)
{
	kh_destroy(arc, p->arc[0]);
	kh_destroy(arc, p->arc[1]);
	free(p->seq); free(p->cov);
}

static void write_node(const fmnode_t *p, long id, kstring_t *out)
{
	int j, first;
	khint_t k;
	if (p->l <= 0) return;
	kputc('@', out); kputl(id, out);
	for (j = 0; j < 2; ++j) {
		kputc('\t', out);
		kputl(p->k[j], out); kputc('>', out);
		if (kh_size(p->arc[j])) {
			for (k = 0, first = 0; k != kh_end(p->arc[j]); ++k) {
				if (kh_exist(p->arc[j], k)) {
					if (first++) kputc(',', out);
					kputl(kh_key(p->arc[j], k), out); kputc(':', out); kputw(kh_val(p->arc[j], k), out);
				}
			}
		} else kputc('.', out);
	}
	kputc('\n', out);
	ks_resize(out, out->l + 2 * p->l + 5);
	for (j = 0; j < p->l; ++j) out->s[out->l++] = "ACGT"[(int)p->seq[j] - 1];
	out->s[out->l] = 0;
	kputsn("\n+\n", 3, out);
	kputsn(p->cov, p->l, out);
	kputc('\n', out);
}

void msg_print(const void *_g)
{
	size_t i;
	const msg_t *g = (const msg_t*)_g;
	kstring_t out;
	out.l = out.m = 0; out.s = 0;
	for (i = 0; i < g->nodes.n; ++i) {
		if (g->nodes.a[i]) {
			out.l = 0;
			if (out.m) out.s[0] = 0;
			write_node(g->nodes.a[i], i, &out);
			fputs(out.s, stdout);
		}
	}
	free(out.s);
}

void *msg_read(const char *fn, int drop_tip)
{
	extern void msg_amend(msg_t *g);
	extern unsigned char seq_nt6_table[128];
	extern void msg_join_unambi(msg_t *g);
	gzFile fp;
	kseq_t *seq;
	khint_t iter;
	size_t i;
	int64_t tot_len = 0;
	double tcpu = cputime();
	msg_t *g;

	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	g = calloc(1, sizeof(msg_t));
	g->h = kh_init(64);
	g->min_ovlp = 1<<30;
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		fmnode_t *p;
		int i, j, ret;
		uint64_t sum, x;
		uint32_t tmp;
		char *q;
		p = msg_node_init();
		p->l = seq->seq.l;
		tmp = p->l + 1;
		kroundup32(tmp);
		p->seq = malloc(tmp); p->cov = malloc(tmp);
		for (i = 0; i < p->l; ++i) p->seq[i] = seq_nt6_table[(int)seq->seq.s[i]];
		strcpy(p->cov, seq->qual.s);
		for (i = 0, sum = 0; i < p->l; ++i) sum += p->cov[i] - 33;
		p->avg_cov = (double)sum / p->l;
		p->aux[0] = p->aux[1] = -1;
		for (j = 0, q = seq->comment.s; j < 2; ++j) {
			p->k[j] = strtol(q, &q, 10);
			if (*q == '>' && q[1] != '.') {
				do {
					++q; x = strtol(q, &q, 10); ++q;
					iter = kh_put(arc, p->arc[j], x, &ret);
					x = strtol(q, &q, 10);
					if (x < g->min_ovlp) g->min_ovlp = x;
					if (ret || kh_val(p->arc[j], iter) < x)
						kh_val(p->arc[j], iter) = x;
				} while (*q == ',');
				++q;
			} else q += 2;
		}
		if (!drop_tip || p->avg_cov > 1.000001 || (kh_size(p->arc[0]) && kh_size(p->arc[1]))) {
			for (j = 0; j < 2; ++j) {
				iter = kh_put(64, g->h, p->k[j], &ret);
				kh_val(g->h, iter) = ret? g->nodes.n<<1|j : (uint64_t)-1;
				if (ret == 0 && fm_verbose >= 2)
					fprintf(stderr, "[W::%s] rank %ld is duplicated.\n", __func__, (long)p->k[j]);
			}
			kv_push(fmnode_p, g->nodes, p);
			tot_len += seq->seq.l;
			if (fm_verbose >= 3 && g->nodes.n % 100000 == 0)
				fprintf(stderr, "[M::%s] read %ld nodes in %ld bp in %.2f sec\n", __func__, g->nodes.n, (long)tot_len, cputime() - tcpu);
		} else msg_node_destroy(p); // an obvious tip or a singleton; do not add
	}
	kseq_destroy(seq);
	gzclose(fp);
	if (fm_verbose >= 2)
		fprintf(stderr, "[M::%s] finished reading %ld nodes in %ld bp in %.2f sec\n", __func__, g->nodes.n, (long)tot_len, cputime() - tcpu);
	msg_amend(g);
	msg_join_unambi(g);
	for (i = 0; i < g->nodes.n; ++i)
		pop_simple_bubble(g, i, 0.1, 2.000001, 1);
	msg_join_unambi(g);
	return g;
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
	int j;
	double tcpu;
	khash_t(arc) *r;
	uint64_t x, z;

	tcpu = cputime(); // excluding time spent on building the hash table
	for (i = 0; i < g->nodes.n; ++i) {
		fmnode_t *p = g->nodes.a[i];
		khint_t k, k2;
		if (p == 0) continue;
		for (j = 0; j < 2; ++j) {
			for (k = 0; k != kh_end(p->arc[j]); ++k) {
				if (!kh_exist(p->arc[j], k)) continue;
				x = kh_key(p->arc[j], k);
				z = get_node_id(g->h, x);
				if (z == (uint64_t)-1) {
					if (fm_verbose >= 5) fprintf(stderr, "[W::%s] rank %ld is non-existing.\n", __func__, (long)x);
					kh_del(arc, p->arc[j], k);
					continue;
				}
				r = g->nodes.a[z>>1]->arc[z&1];
				k2 = kh_get(arc, r, p->k[j]);
				if (k2 == kh_end(r)) {
					if (fm_verbose >= 5) fprintf(stderr, "[W::%s] have %ld->%ld but not the reverse.\n", __func__, (long)p->k[j], (long)x);
					kh_del(arc, p->arc[j], k);
				}
			}
		}
	}
	if (fm_verbose >= 2)
		fprintf(stderr, "[%s] amended the graph in %.2f sec\n", __func__, cputime() - tcpu);
}

static void cut_arc(msg_t *g, uint64_t u, uint64_t v) // delete v from u
{
	uint64_t x;
	khash_t(arc) *q;
	khint_t k;
	x = get_node_id(g->h, u);
	if (x == (uint64_t)-1) return;
	q = g->nodes.a[x>>1]->arc[x&1];
	k = kh_get(arc, q, v);
	if (k != kh_end(q)) kh_del(arc, q, k);
}
/*
static void drop_arc1(fmnode_v *nodes, hash64_t *h, size_t id, int min_ovlp, float min_ovlp_ratio)
{
	fmnode_t *p = &nodes->a[id];
	int j, l, cnt;
	if (min_ovlp == 0 && (min_ovlp_ratio <= 0.01 || min_ovlp_ratio >= 0.99)) return;
	for (j = 0; j < 2; ++j) {
		fm128_v *r = &p->nei[j];
		int max = 0;
		if (r->n == 0) continue;
		for (l = 0; l < r->n; ++l)
			if (r->a[l].y > max) max = r->a[l].y;
		for (l = cnt = 0; l < r->n; ++l) {
			if (r->a[l].y < min_ovlp || (double)r->a[l].y/max < min_ovlp_ratio) {
				if (r->a[l].x != p->k[0] && r->a[l].x != p->k[1])
					cut_arc(nodes, h, r->a[l].x, p->k[j], 1);
				r->a[l].x = 0; // mark the link to delete
				++cnt;
			}
		}
		if (cnt) {
			for (l = cnt = 0; l < r->n; ++l)
				if (r->a[l].x) r->a[cnt++] = r->a[l];
			r->n = cnt;
		}
	}
}
*/

static inline uint64_t first_arc(khash_t(arc) *arc, int *val)
{
	khint_t k;
	for (k = 0; k != kh_end(arc); ++k)
		if (kh_exist(arc, k)) break;
	*val = kh_val(arc, k);
	return kh_key(arc, k);
}

static void add_arc(msg_t *g, uint64_t u, uint64_t v, int ovlp)
{
	uint64_t x = get_node_id(g->h, u);
	khash_t(arc) *q;
	int ret;
	khint_t k;
	if (x == (uint64_t)-1) return;
	q = g->nodes.a[x>>1]->arc[x&1];
	k = kh_put(arc, q, v, &ret);
	if (ret || kh_val(q, k) < ovlp) kh_val(q, k) = ovlp;
}

static void rmnode(msg_t *g, size_t id)
{
	khint_t k0, k1;
	fmnode_t *p = g->nodes.a[id];
	if (p == 0) return;
	if (kh_size(p->arc[0]) && kh_size(p->arc[1])) {
		for (k0 = 0; k0 != kh_end(p->arc[0]); ++k0) {
			if (!kh_exist(p->arc[0], k0)) continue;
			if (kh_key(p->arc[0], k0) == p->k[0] || kh_key(p->arc[0], k0) == p->k[1]) continue; // loop
			for (k1 = 0; k1 != kh_end(p->arc[1]); ++k1) {
				if (!kh_exist(p->arc[1], k1)) continue;
				int ovlp = kh_val(p->arc[0], k0) + kh_val(p->arc[1], k1) - p->l;
				if (ovlp >= g->min_ovlp) {
					add_arc(g, kh_key(p->arc[0], k0), kh_key(p->arc[1], k1), ovlp);
					add_arc(g, kh_key(p->arc[1], k1), kh_key(p->arc[0], k0), ovlp);
				}
			}
		}
	}
	for (k0 = 0; k0 != kh_end(p->arc[0]); ++k0)
		if (kh_exist(p->arc[0], k0)) cut_arc(g, kh_key(p->arc[0], k0), p->k[0]);
	for (k1 = 0; k1 != kh_end(p->arc[1]); ++k1)
		if (kh_exist(p->arc[1], k1)) cut_arc(g, kh_key(p->arc[1], k1), p->k[1]);
	msg_node_destroy(p);
	g->nodes.a[id] = 0;
}

static void rmtip(msg_t *g, int min_len)
{
	size_t i;
	if (min_len < 0) return;
	for (i = 0; i < g->nodes.n; ++i) {
		fmnode_t *p = g->nodes.a[i];
		if (kh_size(p->arc[0]) && kh_size(p->arc[1])) continue; // not a tip
		if (p->l < min_len) rmnode(g, i);
	}
}
/*
static uint64_t pop_bubble(fmnode_v *nodes, hash64_t *h, size_t id_dir, int max_len, fm64_v *stack)
{
	int beg = 0, end = 1, i, j;
	fmnode_t *p = &nodes->a[id_dir>>1], *q;
	fm128_v *r;
	uint64_t x, z, ret = nodes->n;
	if (p->l < 0 || p->nei[id_dir&1].n < 2) return ret;
	p->aux[id_dir&1] = 0;
	stack->n = 0;
	kv_push(uint64_t, *stack, id_dir);
	while (beg < end) {
		for (i = beg; i < end; ++i) {
			x = stack->a[i];
			p = &nodes->a[x>>1];
			if (p->aux[x&1] >= max_len) continue; // stop searching
			r = &p->nei[x&1];
			for (j = 0; j < r->n; ++j) {
				z = get_node_id(h, r->a[j].x);
				if (z == (uint64_t)-1) continue;
				q = &nodes->a[z>>1];
				if (q->l < 0) continue;
				if (q->aux[z&1] >= 0) goto end_db;
				q->aux[z&1] = p->aux[x&1] - r->a[j].y + q->l;
				kv_push(uint64_t, *stack, z^1);
			}
		}
		beg = end; end = stack->n;
	}
end_db:
	if (beg < end && z != id_dir) { // bubble but not a loop
		double min_cov = 1e100;
		int min_i = -1;
		for (i = 0; i < stack->n; ++i) {
			fmnode_t *p = &nodes->a[stack->a[i]>>1];
			if (p->avg_cov < min_cov) min_cov = p->avg_cov, min_i = i;
		}
		ret = stack->a[min_i]>>1;
	}
	for (i = 0; i < stack->n; ++i) // reset aux[]
		nodes->a[stack->a[i]>>1].aux[0] = nodes->a[stack->a[i]>>1].aux[1] = -1;
	return ret;
}
*/
static void pop_simple_bubble(msg_t *g, size_t id, double min_bub_ratio, double min_bub_cov, int max_bubble_diff)
{
	fmnode_t *p = g->nodes.a[id], *q[2], *top_p, *tmp_p;
	khash_t(arc) *r;
	uint64_t nei[2];
	int j, cnt, val;
	double max_cov;
	uint64_t top_id = (uint64_t)-1;
	khint_t k;

	if (p == 0 || kh_size(p->arc[0]) != 1 || kh_size(p->arc[1]) != 1) return;
	for (j = 0; j < 2; ++j) {
		nei[j] = get_node_id(g->h, first_arc(p->arc[j], &val));
		if (nei[j] == (uint64_t)-1) return;
		q[j] = g->nodes.a[nei[j]>>1];
	}
	if (q[0]->l < 0 || q[1]->l < 0) return; // deleted node
	if (kh_size(q[0]->arc[nei[0]&1]) <= 1 || kh_size(q[1]->arc[nei[1]&1]) <= 1) return; // unmerged or inconsistent arc
	r = q[0]->arc[nei[0]&1];
	for (k = 0, max_cov = 0, cnt = 0; k != kh_end(r); ++k) {
		uint64_t t;
		if (!kh_exist(r, k)) continue;
		t = get_node_id(g->h, kh_key(r, k));
		if (t == (uint64_t)-1) continue;
		tmp_p = g->nodes.a[t>>1];
		if (kh_size(tmp_p->arc[0]) != 1 || kh_size(tmp_p->arc[1]) != 1) continue; // skip this node
		if (t&1) flip(tmp_p, g->h); // s.t. tmp_p is on the same strand as p
		if (first_arc(tmp_p->arc[1], &val) != first_arc(p->arc[1], &val)) continue; // not a multi-edge
		if (tmp_p->avg_cov > max_cov) max_cov = tmp_p->avg_cov, top_id = t;
		++cnt;
	}
	if (cnt < 2) return;
	assert(top_id != (uint64_t)-1);
	top_p = g->nodes.a[top_id>>1];
	for (k = 0; k != kh_end(r); ++k) {
		uint64_t t;
		int l, diff, ml = 0, to_del = 0, beg[2], end[2];
		double cov[2];
		if (!kh_exist(r, k)) continue;
		t = get_node_id(g->h, kh_key(r, k));
		if (t == (uint64_t)-1) continue;
		if (top_id == t) continue; // we do not process the node with the highest coverage
		tmp_p = g->nodes.a[t>>1];
		if (kh_size(tmp_p->arc[0]) != 1 || kh_size(tmp_p->arc[1]) != 1) continue; // skip this node
		if (first_arc(tmp_p->arc[1], &val) != first_arc(p->arc[1], &val)) continue; // not a multi-edge
		// the following is really nasty. A banded global alignment would look much cleaner and work better.
		first_arc(top_p->arc[0], &beg[0]); first_arc(top_p->arc[1], &end[0]); end[0] = top_p->l - end[0];
		first_arc(tmp_p->arc[0], &beg[1]); first_arc(tmp_p->arc[1], &end[1]); end[1] = tmp_p->l - end[1];
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
		if (top_p != tmp_p && abs(abs(end[0]-beg[0]) - abs(end[1]-beg[1])) <= max_bubble_diff) {
			if (diff < max_bubble_diff) to_del = 1;
			else if (end[0]-beg[0] - ml < max_bubble_diff || end[1]-beg[1] - ml < max_bubble_diff) to_del = 1;
		}
		if (cov[0] > 0. && cov[1] > 0. && (cov[1] / cov[0] >= min_bub_ratio || cov[1] >= min_bub_cov)) to_del = 0;
		//fprintf(stderr, "%d, %ld, %d, [%f,%f], (%d,%d;%d,%d)\n", to_del, id, diff, cov[0], cov[1], beg[0], end[0], beg[1], end[1]);
		if (to_del) {
			cut_arc(g, tmp_p->k[0], q[0]->k[nei[0]&1]);
			kh_del(arc, r, k);
			cut_arc(g, tmp_p->k[1], q[0]->k[nei[1]&1]);
			cut_arc(g, q[1]->k[nei[1]&1], tmp_p->k[1]);
			msg_node_destroy(tmp_p);
			g->nodes.a[t>>1] = 0;
		}
	}
}

#define __swap(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))

static void flip(fmnode_t *p, hash64_t *h)
{
	extern void seq_reverse(int l, unsigned char *s);
	extern void seq_revcomp6(int l, unsigned char *s);
	khash_t(arc) *t;
	khint_t k;
	seq_revcomp6(p->l, (uint8_t*)p->seq);
	seq_reverse(p->l, (uint8_t*)p->cov);
	__swap(p->k[0], p->k[1]);
	t = p->arc[0]; p->arc[0] = p->arc[1]; p->arc[1] = t;
	k = kh_get(64, h, p->k[0]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
	k = kh_get(64, h, p->k[1]);
	assert(k != kh_end(h));
	kh_val(h, k) ^= 1;
}

static int merge(msg_t *g, size_t w) // merge w's neighbor to the right-end of w
{
	fmnode_t *p = g->nodes.a[w], *q;
	khint_t kp, kq;
	uint32_t tmp;
	int i, j, new_l, val[2];
	uint64_t key[2], qid;
	if (kh_size(p->arc[1]) != 1) return -1; // not exactly one neighbor; cannot be merged
	key[0] = first_arc(p->arc[1], &val[0]);
	kq = kh_get(64, g->h, key[0]);
	assert(kq < kh_end(g->h));
	q = g->nodes.a[qid = kh_val(g->h, kq)>>1];
	if (p == q) return -4; // this is a loop
	if (kh_size(q->arc[kh_val(g->h, kq)&1]) != 1) return -3; // not exactly one neighbor; cannot be merged
	// we can perform a merge
	if (kh_val(g->h, kq)&1) flip(q, g->h); // a "><" bidirectional arc; flip q
	kp = kh_get(64, g->h, p->k[1]); assert(kp != kh_end(g->h)); // get the iterator to p
	kh_del(64, g->h, kp); kh_del(64, g->h, kq); // remove the two ends of the arc in the hash table
	key[1] = first_arc(q->arc[0], &val[1]);
	assert(p->k[1] == key[1] && q->k[0] == key[0]);
	assert(val[0] == val[1]);
	assert(p->l >= val[0] && q->l >= val[0]); // "==" may happen due to end trimming
	new_l = p->l + q->l - val[0];
	tmp = p->l;
	kroundup32(tmp); // this the maximum memory allocated
	if (new_l + 1 >= tmp) { // then double p->seq and p->cov
		tmp = new_l + 1;
		kroundup32(tmp);
		p->seq = realloc(p->seq, tmp);
		p->cov = realloc(p->cov, tmp);
	}
	// merge seq and cov
	for (i = p->l - val[0], j = 0; j < q->l; ++i, ++j) { // write seq and cov
		p->seq[i] = q->seq[j];
		if (i < p->l) {
			if ((int)p->cov[i] + (q->cov[j] - 33) > 126) p->cov[i] = 126;
			else p->cov[i] += q->cov[j] - 33;
		} else p->cov[i] = q->cov[j];
	}
	p->seq[new_l] = p->cov[new_l] = 0;
	p->avg_cov = (p->l * p->avg_cov + q->l * q->avg_cov) / new_l; // recalculate coverage
	p->l = new_l;
	// move neighbors and the rank
	p->k[1] = q->k[1];
	kh_destroy(arc, p->arc[1]);
	p->arc[1] = q->arc[1]; q->arc[1] = 0;
	// update the hash table
	kp = kh_get(64, g->h, p->k[1]);
	assert(kp != kh_end(g->h));
	kh_val(g->h, kp) = w<<1|1;
	// clean up q
	msg_node_destroy(q);
	g->nodes.a[qid] = 0;
	return 0;
}

void msg_join_unambi(msg_t *g)
{
	size_t i;
	double tcpu = cputime();
	for (i = 0; i < g->nodes.n; ++i) {
		if (g->nodes.a[i] == 0) continue;
		while (merge(g, i) == 0);
		flip(g->nodes.a[i], g->h);
		while (merge(g, i) == 0);
	}
	if (fm_verbose >= 2)
		fprintf(stderr, "[%s] joined unambiguous arcs in %.2f sec\n", __func__, cputime() - tcpu);
}
/*
void msg_clean(msg_t *g, const fmclnopt_t *opt)
{
	hash64_t *h = (hash64_t*)g->h;
	fmnode_v *nodes = &g->nodes;
	fm64_v stack;
	size_t i;
	int j;

	kv_init(stack);
	for (j = 0; j < opt->n_iter; ++j) {
		double r = opt->n_iter == 1? 1. : .5 + .5 * j / (opt->n_iter - 1);
		for (i = 0; i < nodes->n; ++i)
			drop_arc1(nodes, h, i, opt->min_ovlp * r, opt->min_ovlp_ratio * r);
		rmtip(nodes, h, opt->min_tip_len * r);
		msg_join_unambi(g);
		for (i = 0; i < nodes->n; ++i) {
			uint64_t x = pop_bubble(nodes, h, i<<1, MAX_POP_EXTENSION, &stack);
			if (x != nodes->n && nodes->a[x].avg_cov < opt->min_weak_cov && nodes->a[x].l < opt->min_tip_len) rmnode(nodes, h, x);
			x = pop_bubble(nodes, h, i<<1|1, MAX_POP_EXTENSION, &stack);
			if (x != nodes->n && nodes->a[x].avg_cov < opt->min_weak_cov && nodes->a[x].l < opt->min_tip_len) rmnode(nodes, h, x);
		}
		msg_join_unambi(g);
		for (i = 0; i < nodes->n; ++i) {
			rm_dup_arc(&nodes->a[i]);
			msg_shrink_node(nodes, h, i);
		}
	}
	for (i = 0; i < nodes->n; ++i)
		if (nodes->a[i].avg_cov < opt->min_weak_cov/2. && nodes->a[i].l < opt->min_tip_len)
			rmnode(nodes, h, i);
	if (opt->min_bub_cov >= 1. && opt->min_bub_ratio < 1.) {
		for (i = 0; i < nodes->n; ++i)
			debubble1_simple(nodes, h, i, opt->min_bub_ratio, opt->min_bub_cov);
		rmtip(nodes, h, opt->min_tip_len);
		msg_join_unambi(g);
	}
	kh_destroy(64, h);
	free(stack.a);
}
*/
