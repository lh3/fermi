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
KHASH_MAP_INIT_INT64(64, fm64_v)

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
			kputl(p->nei[j].a[k].x, out); kputc(':', out);
			kputw((int32_t)p->nei[j].a[k].y, out); kputc(':', out);
			kputw((int32_t)(p->nei[j].a[k].y>>32), out);
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
					x.y = strtol(q, &q, 10); ++q;
					x.y |= (uint64_t)strtol(q, &q, 10)<<32;
					kv_push(fm128_t, p->nei[j], x);
				} while (*q == ',');
				++q;
			} else q += 2;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	return nodes;
}

static inline void put(hash64_t *h, uint64_t key, uint64_t v)
{
	int ret;
	khint_t k;
	k = kh_put(64, h, key, &ret);
	if (ret) kv_init(kh_val(h, k));
	kv_push(uint64_t, kh_val(h, k), v);
}

hash64_t *msg_hash(const fmnode_v *nodes)
{
	size_t i;
	int j, k;
	hash64_t *h;
	h = kh_init(64);
	for (i = 0; i < nodes->n; ++i) {
		const fmnode_t *p = &nodes->a[i];
		if (p->l < 0) continue;
		for (j = 0; j < 2; ++j) {
			put(h, p->k[j], (uint64_t)i<<2 | 0<<1 | j);
			for (k = 0; k < p->nei[j].n; ++k)
				put(h, p->nei[j].a[k].x, (uint64_t)i<<2 | 1<<1 | j);
		}
	}
	return h;
}

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
		if (p->k[0] == 992352 || p->k[1] == 992352) {
			fprintf(stderr, "here\n");
		}
		if (p->nei[0].n && p->nei[1].n) continue; // not a tip
		if (p->avg_cov < min_cov || p->l < min_len) p->l = -1;
		for (j = 0; j < 2; ++j) rmnode(nodes, h, p->k[j]);
	}
}

void msg_rmtip(fmnode_v *nodes, float min_cov, int min_len)
{
	hash64_t *h;
	h = msg_hash(nodes);
	rmtip_core(nodes, h, min_cov, min_len);
}
