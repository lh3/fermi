#include "priv.h"
#include "kvec.h"
#include "kstring.h"

typedef struct {
	int start, len;
	const uint8_t *q;
	const rld_t *e;
	fmintv_v tmpvec[2], match;
} fmsmem_i;

int fm6_smem1_core(const rld_t *e, int len, const uint8_t *q, int x, fmintv_v *mem, int self_match, fmintv_v *prev, fmintv_v *curr)
{
	int i, j, c, ret;
	fmintv_t ik, ok[6];
	fmintv_v *swap;

	fm6_set_intv(e, q[x], ik);

	ik.info = x + 1;
	for (i = x + 1; i < len; ++i) { // forward search
		c = fm6_comp(q[i]);
		fm6_extend(e, &ik, ok, 0);
		if (ok[c].x[2] != ik.x[2]) { // change of the interval size
			if (ik.x[2] != ok[0].x[2]) kv_push(fmintv_t, *curr, ik);
			if (!self_match && ok[0].x[2]) { // some sequences come to an end
				ok[0].info = i;
				kv_push(fmintv_t, *curr, ok[0]);
			}
		}
		if ((!self_match && ok[c].x[2] == 0) || (self_match && ok[c].x[2] < 2)) break; // cannot be extended
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) kv_push(fmintv_t, *curr, ik); // push the last interval if we reach the end
	fm_reverse_fmivec(curr); // s.t. smaller intervals visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;
//	for (i = 0; i < prev->n; ++i) printf("[%lld, %lld, %lld], %lld\n", prev->a[i].x[0], prev->a[i].x[1], prev->a[i].x[2], prev->a[i].info);

	mem->n = 0;
	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			int cont, fl_match; // whether this leads to a full-length read match
			fm6_extend(e, p, ok, 1);
			fl_match = (ok[0].x[2] && p->x[1] < e->mcnt[1]);
			cont = self_match? (ok[c].x[2] > 1) : (ok[c].x[2] != 0);
			if (!cont || fl_match || i == -1) { // keep the hit if: full-length match, reaching the beginning or not extended further
				if (curr->n == 0 || fl_match) { // curr->n to make sure there is no longer matches
//					printf("%d, %lld, [%lld,%lld,%lld]\n", i+1, p->info, p->x[0], p->x[1], p->x[2]);
					if (fl_match || mem->n == 0 || i + 1 < (mem->a[mem->n-1].info>>32&FM_MASK30)) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(ok[0].x[2] != 0) << 63 | (uint64_t)(i + 1)<<32; // bit 64 keeps whether the left-end is closed
						kv_push(fmintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			}
			if (cont && (p->x[1] < e->mcnt[1] || curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = p->info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	fm_reverse_fmivec(mem); // s.t. sorted by the start coordinate

//	for (i = 0; i < mem->n; ++i) printf("[%lld,%lld,%lld] %lld, %lld\n", mem->a[i].x[0], mem->a[i].x[1], mem->a[i].x[2], mem->a[i].info>>32&FM_MASK30, mem->a[i].info&FM_MASK30);
	return ret;
}

fmsmem_i *fm6_miter_init(const rld_t *e, int l, const uint8_t *q)
{
	fmsmem_i *m;
	m = calloc(1, sizeof(fmsmem_i));
	m->len = l, m->e = e, m->q = q;
	return m;
}

void fm6_miter_destroy(fmsmem_i *m)
{
	free(m->tmpvec[0].a); free(m->tmpvec[1].a); free(m->match.a);
	free(m);
}

int fm6_miter_next(fmsmem_i *m)
{
	m->tmpvec[0].n = m->tmpvec[1].n = m->match.n = 0;
	if (m->start >= m->len || m->start < 0) return -1;
	m->start = fm6_smem1_core(m->e, m->len, m->q, m->start, &m->match, 0, &m->tmpvec[0], &m->tmpvec[1]);
	return m->start;
}

int fm6_smem1(const rld_t *e, int len, const uint8_t *q, int x, fmintv_v *mem, int self_match)
{
	int ret;
	fmintv_v a[2];
	kv_init(a[0]); kv_init(a[1]);
	ret = fm6_smem1_core(e, len, q, x, mem, self_match, &a[0], &a[1]);
	free(a[0].a); free(a[1].a);
	return ret;
}

/****************
 ****************/

#include <zlib.h>
#include <pthread.h>
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_DECLARE(64, uint64_t, uint64_t)

typedef khash_t(64) hash64_t;

static uint8_t *paircov(const rld_t *e, int len, const uint8_t *q, int skip, const uint64_t *sorted, hash64_t *h, int *n_supp)
{
	fmsmem_i *iter;
	uint8_t *cov, *pcv;
	cov = calloc((len + 1) * 2, 1);
	pcv = cov + len + 1;
	iter = fm6_miter_init(e, len, q);
	*n_supp = 0;
	while (fm6_miter_next(iter) >= 0) {
		int i, ret, tmp, j;
		uint64_t k, l;
		khint_t kk;
		for (i = 0; i < iter->match.n; ++i) {
			fmintv_t *p = &iter->match.a[i];
			if (p->info>>63 && p->x[1] < e->mcnt[1]) { // full-length match
				tmp = p->info&FM_MASK30;
				for (j = p->info>>32&FM_MASK30; j < tmp; ++j) // update coverage
					if (cov[j] < 255) ++cov[j];
				++(*n_supp);
				for (l = 0; l < p->x[2]; ++l) {
					k = sorted[p->x[1] + l] >> 2;
					if ((k&1) == 0) { // reverse stand; check
						int beg, end;
						kk = kh_get(64, h, k);
						//printf("X %lld get %d\n", k, (int)(p->info&FM_MASK30) - skip);
						if (kk == kh_end(h)) continue; // mate not found on the forward strand
						beg = kh_val(h, kk);
						end = (int)(p->info&FM_MASK30) - skip;
						//printf("%d\t%d\n", beg, end);
						if (beg > end) tmp = beg, beg = end, end = tmp;
						if (end - beg >= FM6_MAX_ISIZE) continue;
						for (j = beg; j < end; ++j)
							if (pcv[j] < 255) ++pcv[j]; // update paired coverage
						kh_del(64, h, kk);
					} else { // forward strand; add
						kk = kh_put(64, h, k^3, &ret);
						tmp = (p->info>>32&FM_MASK30) + skip;
						if (tmp < len) kh_val(h, kk) = tmp;
						else kh_del(64, h, kk);
						//printf("X %lld put %lld, %d<%d\n", k, k^3, tmp, len);
					}
				}
			}
		}
	}
	fm6_miter_destroy(iter);
	kh_clear(64, h);
	return cov;
}

static void paircut_all(const rld_t *e, const uint64_t *sorted, int skip, int n, int *len, char **s, int start, int step, char *const* name)
{
	int i, j;
	hash64_t *h;
	h = kh_init(64);
	for (i = start; i < n; i += step) {
		uint8_t *pcv, *cov, *si = (uint8_t*)s[i];
		int l = len[i], n_supp;
		for (j = 0; j < l; ++j)
			si[j] = seq_nt6_table[si[j]];
		cov = paircov(e, l, si, skip, sorted, h, &n_supp);
		pcv = cov + l + 1;
		for (j = 0; j < l; ++j) {
			si[j] = pcv[j]? "$ACGTN"[si[j]] : "$acgtn"[si[j]];
			cov[j] = cov[j] + 33 < 126? cov[j] + 33 : 126;
		}
		fprintf(stdout, "@%s %d %d\n", name[i], n_supp, l);
		fwrite(si, 1, l, stdout); fwrite("\n+\n", 1, 3, stdout);
		fwrite(cov, 1, l, stdout); fputc('\n', stdout);
		free(cov);
	}
	kh_destroy(64, h);
}

typedef struct {
	int n, m, *l;
	char **s, **name;
} seqbuf_t;

static int fill_seqbuf(kseq_t *kseq, seqbuf_t *buf, int64_t max_len)
{
	int64_t l = 0;
	int i;
	for (i = 0; i < buf->n; ++i) {
		free(buf->s[i]);
		free(buf->name[i]);
	}
	buf->n = 0;
	while (kseq_read(kseq) >= 0) {
		if (buf->n == buf->m) {
			buf->m = buf->m? buf->m<<1 : 256;
			buf->l = realloc(buf->l, buf->m * sizeof(int));
			buf->s = realloc(buf->s, buf->m * sizeof(void*));
			buf->name = realloc(buf->name, buf->m * sizeof(void*));
		}
		buf->l[buf->n] = kseq->seq.l;
		buf->s[buf->n] = strdup(kseq->seq.s);
		buf->name[buf->n++] = strdup(kseq->name.s);
		l += kseq->seq.l;
		if (l >= max_len) break;
	}
	return buf->n;
}

typedef struct {
	const rld_t *e;
	const uint64_t *sorted;
	int start, step, skip;
	seqbuf_t *buf;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	paircut_all(w->e, w->sorted, w->skip, w->buf->n, w->buf->l, w->buf->s, w->start, w->step, w->buf->name);
	return 0;
}

int fm6_paircut(const char *fn, const rld_t *e, uint64_t *sorted, int skip, int n_threads)
{
	int i;
	kseq_t *seq;
	gzFile fp;
	worker_t *w;
	seqbuf_t *buf;
	pthread_t *tid;
	pthread_attr_t attr;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	buf = calloc(1, sizeof(seqbuf_t));
	for (i = 0; i < n_threads; ++i)
		w[i].e = e, w[i].sorted = sorted, w[i].skip = skip, w[i].step = n_threads, w[i].start = i, w[i].buf = buf;
	
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);

	while (fill_seqbuf(seq, buf, 1<<28) > 0) {
		for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], &attr, worker, w + i);
		for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	}

	free(buf->l); free(buf->s); free(buf->name);
	kseq_destroy(seq);
	gzclose(fp);
	free(tid); free(w);
	return 0;
}

// legacy interface
int fm6_smem(const rld_t *e, int len, const uint8_t *q, fmintv_v *mem, int self_match)
{
	int x = 0, i;
	fmintv_v tmp;
	kv_init(tmp);
	mem->n = 0;
	do {
		x = fm6_smem1(e, len, q, x, &tmp, self_match);
		for (i = 0; i < tmp.n; ++i) {
			kv_push(fmintv_t, *mem, tmp.a[i]);
		}
	} while (x < len);
	return mem->n;
}

int fm6_write_smem(const rld_t *e, const fmintv_t *a, kstring_t *s)
{
	s->l = 0;
	kputuw(a->info>>32&FM_MASK30, s); kputc('\t', s); kputuw(a->info&FM_MASK30, s); kputc('\t', s);
	kputuw(a->x[2] > 0xffffffffU? 0xffffffffU : a->x[2], s); kputc('\t', s);
	kputc("OT"[a->info>>63], s); kputc("OT"[a->x[1] < e->mcnt[1]], s);
	return s->l;
}
