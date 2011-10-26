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

typedef khash_t(64) hash64_t;

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

static void check(fmnode_v *nodes, hash64_t *h)
{
	size_t i;
	int j, l, ll;
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
					continue;
				}
				q = &nodes->a[kh_val(h, k)>>1];
				r = &q->nei[kh_val(h, k)&1];
				for (ll = 0; ll < r->n; ++ll)
					if (r->a[ll].x == p->k[j]) break;
				if (ll == r->n) {
					if (fm_verbose >= 2) fprintf(stderr, "[W::%s] have %ld->%ld but no reverse.\n", __func__, (long)p->k[j], (long)x);
					continue;
				}
			}
		}
	}
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

static void rmtip(fmnode_v *nodes, hash64_t *h, float min_cov, int min_len)
{
	size_t i;
	int j, l;
	for (i = 0; i < nodes->n; ++i) {
		fmnode_t *p = &nodes->a[i];
		if (p->nei[0].n && p->nei[1].n) continue; // not a tip
		if (p->avg_cov < min_cov && p->l < min_len) {
			p->l = -1;
			for (j = 0; j < 2; ++j)
				for (l = 0; l < p->nei[j].n; ++l)
					cut_arc(nodes, h, p->nei[j].a[l].x, p->k[j], 1);
		}
	}
}

static void debubble1_simple(fmnode_v *nodes, hash64_t *h, size_t id, hash64_t *tmph)
{
	fmnode_t *p = &nodes->a[id], *q[2];
	fm128_v *r[2];
	uint64_t nei[2];
	int j, ret;
	khint_t k;
	if (p->l <= 0 || p->nei[0].n != 1 || p->nei[1].n != 1) return;
	for (j = 0; j < 2; ++j) {
		nei[j] = get_node_id(h, p->nei[j].a[0].x);
		q[j] = &nodes->a[nei[j]>>1];
	}
	if (q[0]->l < 0 || q[1]->l < 0) return; // deleted node
	if (q[0]->nei[nei[0]&1].n <= 1 || q[1]->nei[nei[1]&1].n <= 1) return; // unmerged or inconsistent arc
	if (q[0]->nei[nei[0]&1].n != q[1]->nei[nei[1]&1].n) return; // branching
	r[0] = &q[0]->nei[nei[0]&1]; r[1] = &q[1]->nei[nei[1]&1];
	kh_clear(64, tmph);
	for (j = 0; j < r[0]->n; ++j) {
		k = kh_put(64, tmph, r[0]->a[j].x, &ret);
		assert(ret == 0);
	}
	for (j = 0; j < r[1]->n; ++j) {
		k = kh_get(64, tmph, r[1]->a[j].x);
		if (k == kh_end(tmph)) break;
	}
	if (j != r[1]->n) return; // branching
}

static void erode_end1(fmnode_v *nodes, hash64_t *h, size_t id, int min_cov)
{
	fmnode_t *p = &nodes->a[id];
	int j, l, k, el[2], ol[2];
	if (p->l <= 0) return;
	for (j = 0; j < 2; ++j) {
		for (l = 0, ol[j] = 0; l < p->nei[j].n; ++l)
			if (ol[j] < p->nei[j].a[l].y) ol[j] = p->nei[j].a[l].y;
	}
	for (l = 0; l < p->l - ol[1]; ++l)
		if (p->cov[l] - 33 > min_cov) break;
	el[0] = l;
	for (l = p->l - 1; l >= ol[0]; --l)
		if (p->cov[l] - 33 > min_cov) break;
	el[1] = p->l - 1 - l;
	if (el[0] + el[1] >= p->l) return; // erode the entire node
	for (j = 0; j < 2; ++j) {
		if (el[j] == 0) continue;
		for (l = 0; l < p->nei[j].n; ++l) {
			fm128_t *q = &p->nei[j].a[l];
			fmnode_t *r;
			uint64_t nei;
			nei = get_node_id(h, q->x);
			if (nei == (uint64_t)-1) continue;
			q->y = q->y < el[j]? 0 : q->y - el[j];
			r = &nodes->a[nei>>1];
			for (k = 0; k < r->nei[nei&1].n; ++k)
				if (r->nei[nei&1].a[k].x == p->k[j])
					break;
			if (k < r->nei[nei&1].n) // no broken edges
				r->nei[nei&1].a[k].y = q->y;
		}
	}
	if (el[1]) p->l -= el[1];
	if (el[0]) {
		p->l -= el[0];
		memmove(p->cov, p->cov + el[0], p->l);
		memmove(p->seq, p->seq + el[0], p->l);
	}
	p->cov[p->l] = p->seq[p->l] = 0;
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
	hash64_t *h, *tmph;
	size_t i;
	tmph = kh_init(64);
	h = build_hash(nodes);
	if (opt->check) check(nodes, h);
	if (opt->min_tip_len && opt->min_tip_cov >= 1.)
		for (i = 0; i < 5; ++i)
			rmtip(nodes, h, opt->min_tip_cov, opt->min_tip_len);
	merge(nodes, h, 0);
	if (opt->min_term_cov)
		for (i = 0; i < nodes->n; ++i)
			erode_end1(nodes, h, i, opt->min_term_cov);
	clean_core(nodes, h);
	kh_destroy(64, h);
	kh_destroy(64, tmph);
}
