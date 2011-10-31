#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fermi.h"
#include "rld.h"
#include "kstring.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khash.h"
KHASH_MAP_INIT_INT64(64, uint64_t)

#define fm128_xlt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y > (b).y))
#include "ksort.h"
KSORT_INIT(128x, fm128_t, fm128_xlt)

#define MAX_DEBUBBLE_DIFF 3
#define MAX_POP_EXTENSION 1000

typedef khash_t(64) hash64_t;

static void flip(fmnode_t *p, hash64_t *h);

void msg_write_node(const fmnode_t *p, long id, kstring_t *out)
{
	int j, k;
	if (p->l <= 0) return;
	kputc('@', out); kputl(id, out);
	for (j = 0; j < 2; ++j) {
		kputc('\t', out);
		kputl(p->k[j], out); kputc('>', out);
		for (k = 0; k < p->nei[j].n; ++k) {
			if (k) kputc(',', out);
			kputl(p->nei[j].a[k].x, out); kputc(':', out); kputw((int32_t)p->nei[j].a[k].y, out);
			//kputc(':', out); kputw((int32_t)(p->nei[j].a[k].y>>32), out);
		}
		if (p->nei[j].n == 0) kputc('.', out);
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

fmnode_v *msg_read(const char *fn)
{
	extern unsigned char seq_nt6_table[128];
	gzFile fp;
	kseq_t *seq;
	fmnode_v *nodes;
	int64_t cnt = 0, tot_len = 0;

	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	seq = kseq_init(fp);
	nodes = calloc(1, sizeof(fmnode_v));
	while (kseq_read(seq) >= 0) {
		fmnode_t *p;
		int i, j;
		uint64_t sum;
		uint32_t tmp;
		char *q;
		kv_pushp(fmnode_t, *nodes, &p);
		kv_init(p->nei[0]); kv_init(p->nei[1]);
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
			p->k[j] = strtol(q, &q, 10);
			if (*q == '>' && q[1] != '.') {
				fm128_t x;
				do {
					++q;
					x.x = strtol(q, &q, 10); ++q;
					x.y = strtol(q, &q, 10);
					kv_push(fm128_t, p->nei[j], x);
				} while (*q == ',');
				++q;
			} else q += 2;
		}
		msg_rmdup(p);
		++cnt; tot_len += seq->seq.l;
		if (fm_verbose >= 3 && cnt % 100000 == 0)
			fprintf(stderr, "[%s] read %ld nodes in %ld bp\n", __func__, (long)cnt, (long)tot_len);
	}
	kseq_destroy(seq);
	gzclose(fp);
	if (fm_verbose >= 3)
		fprintf(stderr, "[%s] In total: %ld nodes in %ld bp\n", __func__, (long)cnt, (long)tot_len);
	return nodes;
}

void msg_rmdup(fmnode_t *p)
{
	int j, l, cnt;
	uint64_t x;
	for (j = 0; j < 2; ++j) {
		fm128_v *r = &p->nei[j];
		if (r->n < 2) continue;
		ks_introsort(128x, r->n, r->a);
		x = r->a[0].x;
		for (l = 1, cnt = 0; l < r->n; ++l) {
			if (r->a[l].x == x) r->a[l].x = 0, ++cnt;
			else x = r->a[l].x;
		}
		if (cnt) {
			for (l = 0, cnt = 0; l < r->n; ++l)
				if (r->a[l].x) r->a[cnt++] = r->a[l];
			r->n = cnt;
		}
	}
}

static hash64_t *build_hash(const fmnode_v *nodes)
{
	size_t i;
	int j, ret;
	hash64_t *h;
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
	return h;
}

static inline uint64_t get_node_id(hash64_t *h, uint64_t tid)
{
	khint_t k;
	k = kh_get(64, h, tid);
	if (k == kh_end(h)) return (uint64_t)(-1);
	return kh_val(h, k);
}

static int check(fmnode_v *nodes, hash64_t *h)
{
	size_t i;
	int j, l, ll, ret = 0;
	khint_t k;
	for (i = 0; i < nodes->n; ++i) {
		fmnode_t *q, *p = &nodes->a[i];
		fm128_v *r;
		if (p->l <= 0) continue;
		for (j = 0; j < 2; ++j) {
			for (l = 0; l < p->nei[j].n; ++l) {
				uint64_t x = p->nei[j].a[l].x;
				k = kh_get(64, h, x);
				if (k == kh_end(h)) {
					if (fm_verbose >= 2) fprintf(stderr, "[W::%s] tip %ld is non-existing.\n", __func__, (long)x);
					ret = -1;
					continue;
				}
				q = &nodes->a[kh_val(h, k)>>1];
				r = &q->nei[kh_val(h, k)&1];
				for (ll = 0; ll < r->n; ++ll)
					if (r->a[ll].x == p->k[j]) break;
				if (ll == r->n) {
					if (fm_verbose >= 2) fprintf(stderr, "[W::%s] have %ld->%ld but no reverse.\n", __func__, (long)p->k[j], (long)x);
					return -1;
					continue;
				}
			}
		}
	}
	return ret;
}

static void cut_arc(fmnode_v *nodes, hash64_t *h, uint64_t u, uint64_t v, int remove) // delete v from u
{
	khint_t k;
	int i, j;
	fmnode_t *p;
	fm128_v *r;
	k = kh_get(64, h, u);
	if (k == kh_end(h)) return;
	p = &nodes->a[kh_val(h, k)>>1];
	r = &p->nei[kh_val(h, k)&1];
	if (remove) {
		for (j = i = 0; j < r->n; ++j)
			if (r->a[j].x != v) r->a[i++] = r->a[j];
		r->n = i;
	} else {
		for (j = 0; j < r->n; ++j)
			if (r->a[j].x == v) r->a[j].x = r->a[j].y = 0;
	}
}

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
				if (r->a[l].x != p->k[j])
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

static inline void rmnode(fmnode_v *nodes, hash64_t *h, size_t id)
{
	int j, l;
	fmnode_t *p = &nodes->a[id];
	p->l = -1;
	for (j = 0; j < 2; ++j)
		for (l = 0; l < p->nei[j].n; ++l)
			cut_arc(nodes, h, p->nei[j].a[l].x, p->k[j], 1);
}

static void rmtip(fmnode_v *nodes, hash64_t *h, int min_len)
{
	size_t i;
	if (min_len < 0) return;
	for (i = 0; i < nodes->n; ++i) {
		fmnode_t *p = &nodes->a[i];
		if (p->nei[0].n && p->nei[1].n) continue; // not a tip
		if (p->l < min_len) rmnode(nodes, h, i);
	}
}

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

static void debubble1_simple(fmnode_v *nodes, hash64_t *h, size_t id, double min_bub_ratio, double min_bub_cov)
{
	fmnode_t *p = &nodes->a[id], *q[2], *top_p, *tmp_p;
	fm128_v *r[2];
	uint64_t nei[2];
	int j, cnt;
	double max_cov;
	uint64_t top_id = (uint64_t)-1;

	if (p->l <= 0 || p->nei[0].n != 1 || p->nei[1].n != 1) return;
	for (j = 0; j < 2; ++j) {
		nei[j] = get_node_id(h, p->nei[j].a[0].x);
		q[j] = &nodes->a[nei[j]>>1];
	}
	if (q[0]->l < 0 || q[1]->l < 0) return; // deleted node
	if (q[0]->nei[nei[0]&1].n <= 1 || q[1]->nei[nei[1]&1].n <= 1) return; // unmerged or inconsistent arc
	r[0] = &q[0]->nei[nei[0]&1]; r[1] = &q[1]->nei[nei[1]&1];
	for (j = 0, max_cov = 0, cnt = 0; j < r[0]->n; ++j) {
		uint64_t t = get_node_id(h, r[0]->a[j].x);
		tmp_p = &nodes->a[t>>1];
		if (tmp_p->nei[0].n != 1 || tmp_p->nei[1].n != 1) continue; // skip this node
		if (t&1) flip(tmp_p, h); // s.t. tmp_p is on the same strand as p
		if (tmp_p->nei[1].a[0].x != p->nei[1].a[0].x) continue; // not a multi-edge
		if (tmp_p->avg_cov > max_cov) max_cov = tmp_p->avg_cov, top_id = t;
		++cnt;
	}
	if (cnt < 2) return;
	assert(top_id != (uint64_t)-1);
	top_p = &nodes->a[top_id>>1];
	for (j = 0; j < r[0]->n; ++j) {
		uint64_t t = get_node_id(h, r[0]->a[j].x);
		int l, diff, ml = 0, to_del = 0, beg[2], end[2];
		double cov[2];
		if (top_id == t) continue; // we do not process the node with the highest coverage
		tmp_p = &nodes->a[t>>1];
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
			cut_arc(nodes, h, tmp_p->k[0], q[0]->k[nei[0]&1], 1);
			cut_arc(nodes, h, q[0]->k[nei[0]&1], tmp_p->k[0], 0);
			cut_arc(nodes, h, tmp_p->k[1], q[0]->k[nei[1]&1], 1);
			cut_arc(nodes, h, q[1]->k[nei[1]&1], tmp_p->k[1], 1);
			tmp_p->l = -1;
		}
	}
	for (j = 0, cnt = 0; j < r[0]->n; ++j)
		if (r[0]->a[j].x) r[0]->a[cnt++] = r[0]->a[j];
	r[0]->n = cnt;
}

#define __swap(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))

static void flip(fmnode_t *p, hash64_t *h)
{
	extern void seq_reverse(int l, unsigned char *s);
	extern void seq_revcomp6(int l, unsigned char *s);
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
}

static int merge(fmnode_v *nodes, hash64_t *h, size_t w) // merge i's neighbor to the right-end of i
{
	fmnode_t *p = &nodes->a[w], *q;
	khint_t kp, kq;
	uint32_t tmp;
	int i, j, new_l;
	if (p->nei[1].n != 1) return -1; // cannot be merged
	kq = kh_get(64, h, p->nei[1].a[0].x);
	if (kq == kh_end(h)) return -2; // not found
	q = &nodes->a[kh_val(h, kq)>>1];
	if (p == q) return -4; // this is a loop
	if (q->nei[kh_val(h, kq)&1].n != 1) return -3; // cannot be merged
	// we can perform a merge
	if (kh_val(h, kq)&1) flip(q, h); // a "><" bidirectional arc; flip q
	kp = kh_get(64, h, p->k[1]); assert(kp != kh_end(h)); // get the iterator to p
	kh_del(64, h, kp); kh_del(64, h, kq); // remove the two ends of the arc in the hash table
	assert(p->k[1] == q->nei[0].a[0].x);
	assert(q->k[0] == p->nei[1].a[0].x);
	assert(p->nei[1].a[0].y == q->nei[0].a[0].y);
	assert(p->l >= p->nei[1].a[0].y && q->l >= p->nei[1].a[0].y); // "==" may happen due to end trimming
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
	// update the hash table
	kp = kh_get(64, h, p->k[1]);
	assert(kp != kh_end(h));
	kh_val(h, kp) = w<<1|1;
	// clean up q
	free(q->cov); free(q->seq);
	q->cov = q->seq = 0; q->l = -1;
	free(q->nei[0].a);
	q->nei[0].n = q->nei[0].m = q->nei[1].n = q->nei[1].m = 0;
	q->nei[0].a = q->nei[1].a = 0; // q->nei[1] has been copied over to p->nei[1], so we can delete it
	return 0;
}

static void clean_core(fmnode_v *nodes, hash64_t *h)
{
	size_t i;
	for (i = 0; i < nodes->n; ++i) {
		if (nodes->a[i].l <= 0) continue;
		while (merge(nodes, h, i) == 0);
		flip(&nodes->a[i], h);
		while (merge(nodes, h, i) == 0);
	}
}

void msg_clean(fmnode_v *nodes, const fmclnopt_t *opt)
{
	hash64_t *h;
	fm64_v stack;
	size_t i;
	int j;
	kv_init(stack);
	h = build_hash(nodes);
	if (opt->check) check(nodes, h);
	for (j = 0; j < opt->n_iter; ++j) {
		double r = opt->n_iter == 1? 1. : .5 + .5 * j / (opt->n_iter - 1);
		for (i = 0; i < nodes->n; ++i) msg_rmdup(&nodes->a[i]);
		for (i = 0; i < nodes->n; ++i)
			drop_arc1(nodes, h, i, opt->min_ovlp * r, opt->min_ovlp_ratio * r);
		rmtip(nodes, h, opt->min_tip_len * r);
		clean_core(nodes, h);
		for (i = 0; i < nodes->n; ++i) {
			uint64_t x = pop_bubble(nodes, h, i<<1, MAX_POP_EXTENSION, &stack);
			if (x != nodes->n && nodes->a[x].avg_cov < opt->min_weak_cov) rmnode(nodes, h, x);
			x = pop_bubble(nodes, h, i<<1|1, MAX_POP_EXTENSION, &stack);
			if (x != nodes->n && nodes->a[x].avg_cov < opt->min_weak_cov) rmnode(nodes, h, x);
		}
		for (i = 0; i < nodes->n; ++i) msg_rmdup(&nodes->a[i]);
		clean_core(nodes, h);
	}
//	for (i = 0; i < nodes->n; ++i)
//		if (nodes->a[i].avg_cov < 1.01 && nodes->a[i].nei[0].n == 1&& nodes->a[i].nei[1].n) rmnode(nodes, h, i);
//	clean_core(nodes, h);
	if (opt->min_bub_cov >= 1. && opt->min_bub_ratio < 1.) {
		for (i = 0; i < nodes->n; ++i)
			debubble1_simple(nodes, h, i, opt->min_bub_ratio, opt->min_bub_cov);
		rmtip(nodes, h, opt->min_tip_len);
		clean_core(nodes, h);
	}
	kh_destroy(64, h);
	free(stack.a);
}
