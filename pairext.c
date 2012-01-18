#include <zlib.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include "priv.h"
#include "kstring.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#define BUF_SIZE   0x10000
#define MIN_INSERT 50

extern unsigned char seq_nt6_table[128];

typedef struct {
	int len;
	uint8_t *semitig;
	fm64_v reads;
} ext1_t;

static int read_unitigs(kseq_t *kseq, int n, ext1_t *buf, int min_dist, int max_dist)
{
	int k, j = 0;
	fm128_v reads;
	assert(n >= 3);
	kv_init(reads);
	while (kseq_read(kseq) >= 0) {
		char *q;
		ext1_t *p;
		int beg, end;
		if (kseq->comment.l == 0) continue; // no comments
		if (kseq->seq.l < min_dist) continue; // too short; skip
		reads.n = 0;
		for (k = 0, q = kseq->comment.s; *q && k < 2; ++q) // skip the first two fields
			if (isspace(*q)) ++k;
		while (isdigit(*q)) { // read mapping
			fm128_t x;
			x.x = strtol(q, &q, 10); ++q;
			x.y = *q == '-'? 1 : 0; q += 2;
			x.y |= (uint64_t)strtol(q, &q, 10)<<32; ++q;
			x.y |= strtol(q, &q, 10)<<1;
			kv_push(fm128_t, reads, x);
			if (*q++ == 0) break;
		}
		// left-end
		p = &buf[j++];
		p->reads.n = 0;
		end = kseq->seq.l < max_dist? kseq->seq.l : max_dist;
		for (k = 0; k < reads.n; ++k) {
			fm128_t *r = &reads.a[k];
			if ((r->y&1) && r->y<<32>>33 <= end)
				kv_push(uint64_t, p->reads, (r->x^1)<<1);
		}
		if (p->reads.n) { // potentially extensible
			p->len = end;
			p->semitig = calloc(p->len + 1, 1);
			for (k = 0; k < end; ++k)
				p->semitig[k] = seq_nt6_table[(int)kseq->seq.s[k]];
		} else --j;
		// right-end
		p = &buf[j++];
		p->reads.n = 0;
		beg = kseq->seq.l < max_dist? 0 : kseq->seq.l - max_dist;
		for (k = 0; k < reads.n; ++k) {
			fm128_t *r = &reads.a[k];
			if ((r->y&1) == 0 && r->y>>32 >= beg)
				kv_push(uint64_t, p->reads, (r->x^1)<<1|1);
		}
		if (p->reads.n) {
			p->len = kseq->seq.l - beg;
			p->semitig = calloc(p->len + 1, 1);
			for (k = 0; k < p->len; ++k)
				p->semitig[k] = seq_nt6_table[(int)kseq->seq.s[k + beg]];
		} else --j;
		if (j + 1 >= n) break;
	}
	return j;
}

static void pext_core(const rld_t *e, int n, ext1_t *buf, int start, int step, int is_aggressive)
{
	extern void seq_reverse(int l, unsigned char *s);
	kstring_t rd, seq, out;
	int i, j;
	rd.l = seq.l = rd.m = seq.m = out.l = out.m = 0; rd.s = seq.s = out.s = 0;
	for (i = start; i < n; i += step) {
		ext1_t *p = &buf[i];
		msg_t *g;
		int tmp, max_len = 0;
		seq.l = 0;
		kputsn((char*)p->semitig, p->len + 1, &seq); // +1 to include the ending NULL
		for (j = 0; j < p->reads.n; ++j) {
			assert(p->reads.a[j] < e->mcnt[1]);
			fm_retrieve(e, p->reads.a[j], &rd);
			if (rd.l > max_len) max_len = rd.l;
			seq_reverse(rd.l, (uint8_t*)rd.s);
			kputsn(rd.s, rd.l + 1, &seq); // +1 to include NULL
		}
		if (max_len >= p->len) continue; // a read is longer than the semitig? stop; FIXME: this can be improved
		// de novo assembly
		fm6_api_correct(16, seq.l, seq.s, 0);
		g = fm6_api_unitig(-1, seq.l, seq.s);
		msg_rm_tips(g, max_len + 1, 2); // very mild tip removal; note that max_len+1 <= p->len
		msg_join_unambi(g);
		if (is_aggressive) { // aggressive bubble popping
			msg_popbub_open(g, max_len - 1, 1);
			msg_join_unambi(g);
			msg_popbub(g, max_len * 2, 25);
			msg_join_unambi(g);
		}
		// decide if keep the longest contig
		for (j = 0, max_len = 0, tmp = -1; j < g->nodes.n; ++j)
			if (g->nodes.a[j].l >= max_len)
				max_len = g->nodes.a[j].l, tmp = j;
		if (max_len > p->len) { // the semitig is extensible
			fmnode_t *q = &g->nodes.a[tmp];
			out.l = 0;
			kputc('>', &out); kputw(i, &out); kputc(' ', &out); kputw(max_len - p->len, &out); kputc('\n', &out);
			for (j = 0; j < q->l; ++j)
				kputc("$ACGTN"[(int)q->seq[j]], &out);
			puts(out.s);
		}
		msg_destroy(g);
	}
	free(rd.s); free(seq.s); free(out.s);
}

typedef struct {
	const rld_t *e;
	ext1_t *buf;
	int n, start, step, is_aggressive;
} worker_t;

static void *worker(void *_w)
{
	worker_t *w = (worker_t*)_w;
	pext_core(w->e, w->n, w->buf, w->start, w->step, w->is_aggressive);
	return 0;
}

int fm6_pairext(const rld_t *e, const char *fng, int n_threads, double avg, double std, int is_aggressive)
{
	pthread_t *tid;
	pthread_attr_t attr;
	worker_t *w;
	int min_dist, max_dist;
	kseq_t *kseq;
	gzFile fp;
	ext1_t *buf;
	int i, n, old_verbose;

	max_dist = (int)(avg + std * 2. + .499);
	min_dist = (int)(avg - std * 2. + .499);
	if (min_dist < MIN_INSERT) min_dist = MIN_INSERT;
	fp = strcmp(fng, "-")? gzopen(fng, "rb") : gzdopen(fileno(stdin), "rb");
	if (fp == 0) return -1;
	kseq = kseq_init(fp);
	buf = calloc(BUF_SIZE, sizeof(ext1_t));

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) {
		w[i].e = e;
		w[i].buf = buf;
		w[i].start = i;
		w[i].step = n_threads;
		w[i].is_aggressive = is_aggressive;
	}
	old_verbose = fm_verbose; fm_verbose = 1; // to suppress messages and warnings
	while ((n = read_unitigs(kseq, BUF_SIZE, buf, min_dist, max_dist)) > 0) {
		for (i = 0; i < n_threads; ++i) w[i].n = n;
		for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], &attr, worker, w + i);
		for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	}
	fm_verbose = old_verbose;

	for (i = 0; i < BUF_SIZE; ++i) {
		free(buf[i].semitig);
		free(buf[i].reads.a);
	}
	free(buf);
	free(tid); free(w);
	kseq_destroy(kseq);
	gzclose(fp);
	return 0;
}
