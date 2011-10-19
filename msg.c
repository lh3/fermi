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
		char *q;
		kv_pushp(fmnode_t, *nodes, &p);
		kv_init(p->nei[0]); kv_init(p->nei[1]);
		p->l = seq->seq.l;
		p->seq = malloc(p->l);
		for (i = 0; i < p->l; ++i) p->seq[i] = seq_nt6_table[(int)seq->seq.s[i]];
		p->cov = strdup(seq->qual.s);
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

static void amend(fmnode_v *nodes, hash64_t *h)
{
	size_t i;
	int j, l;
	khint_t k;
	for (i = 0; i < nodes->n; ++i) {
		fmnode_t *p = &nodes->a[i];
		for (j = 0; j < 2; ++j) {
			for (l = 0; l < p->nei[j].n; ++l) {
				uint64_t x = p->nei[j].a[l].x;
				k = kh_get(64, h, x);
				if (k == kh_end(h)) {
					if (fm_verbose >= 2) fprintf(stderr, "[W::%s] tip %ld is non-existing.\n", __func__, (long)x);
					continue;
				}
			}
		}
	}
}

/*
static int rmnode(fmnode_v *nodes, hash64_t *h, uint64_t k0)
{
	khint_t k;
	fm64_v *q;
	int i, j, n, cnt = 0;
	k = kh_get(64, h, k0);
	q = &kh_val(h, k);
	for (i = 0; i < q->n; ++i) {
		fm128_v *r;
		if ((q->a[i] & 2) == 0) continue; // not in the nei array
		r = &nodes->a[q->a[i]>>2].nei[q->a[i]&1];
		for (j = n = 0; j < r->n; ++j)
			if (r->a[j].x != k0) r->a[n++] = r->a[j];
		cnt += r->n - n;
		r->n = n;
	}
	kh_del(64, h, k);
	return cnt;
}

static void rmtip_core(fmnode_v *nodes, hash64_t *h, float min_cov, int min_len)
{
	size_t i;
	int j;
	for (i = 0; i < nodes->n; ++i) {
		fmnode_t *p = &nodes->a[i];
		if (p->nei[0].n && p->nei[1].n) continue; // not a tip
		if (p->avg_cov < min_cov || p->l < min_len) p->l = -1;
		for (j = 0; j < 2; ++j) rmnode(nodes, h, p->k[j]);
	}
}

#define __swap(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))

static void flip(fmnode_t *p)
{
	extern void seq_reverse(int l, unsigned char *s);
	extern void seq_revcomp6(int l, unsigned char *s);
	fm128_v t;
	seq_revcomp6(p->l, (uint8_t*)p->seq);
	seq_reverse(p->l, (uint8_t*)p->cov);
	__swap(p->k[0], p->k[1]);
	t = p->nei[0]; p->nei[0] = p->nei[1]; p->nei[1] = t;
}

static int merge_node(fmnode_v *nodes, hash64_t *h, uint64_t key0, uint64_t key1, int len)
{
	int j;
	fmnode_t *p[2];
	fm64_v *q[2];
	khint_t k[2];
	uint32_t key[2];
	key[0] = key0; key[1] = key1;
	for (j = 0; j < 2; ++j) {
		k[j] = kh_get(64, h, key[j]);
		if (k[j] == kh_end(h)) return -1; // cannot find in the hash table
		q[j] = &kh_val(h, k[j]);
		if (q[j]->n != 2) return -1; // not exactly two nodes
		if (((q[j]->a[0] ^ q[j]->a[1]) & 2) == 0) return -1; // not a start and an end
		if (q[j]->a[0]&2) __swap(q[j]->a[0], q[j]->a[1]);
		p[j] = &nodes->a[q[j]->a[0]>>2];
	}
	if (q[0]->a[0]>>2 != q[1]->a[1]>>2 || q[0]->a[1]>>2 != q[1]->a[0]>>2) return -1;
	if (!(q[0]->a[0]&1)) flip(p[0]);
	if (q[0]->a[1]&1) flip(p[1]);
	return 0;
}

static void clean_core(fmnode_v *nodes, hash64_t *h)
{
	size_t i;
	for (i = 0; i < nodes->n; ++i) {
		fmnode_t *p = &nodes->a[i];
		if (p->nei[0].n == 1) merge_node(nodes, h, p->k[0], p->nei[0].a[0].x, (uint32_t)p->nei[0].a[0].y);
		if (p->nei[1].n == 1) merge_node(nodes, h, p->k[1], p->nei[1].a[0].x, (uint32_t)p->nei[1].a[0].y);
	}
}
*/
void msg_clean(fmnode_v *nodes, float min_cov, int min_len)
{
	hash64_t *h;
	h = build_hash(nodes);
	amend(nodes, h);
//	rmtip_core(nodes, h, min_cov, min_len);
//	clean_core(nodes, h);
}
