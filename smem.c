#include <math.h>
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
	if (i == len) { // push the last interval if we reach the end
		kv_push(fmintv_t, *curr, ik); // always push this interval
		if (!self_match) { // then test if this last interval is terminated
			fm6_extend(e, &ik, ok, 0);
			if (ok[0].x[2]) { // terminated; then push
				ok[0].info = len;
				kv_push(fmintv_t, *curr, ok[0]);
			}
		}
	}
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

/***************
 * Remap reads *
 ***************/

#include <zlib.h>
#include <ctype.h>
#include <pthread.h>
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#include "khash.h"
KHASH_DECLARE(64, uint64_t, uint64_t)

typedef khash_t(64) hash64_t;

typedef struct {
	int n_supp, len;
	uint8_t *cov, *pcv; // cov and pcv are allocated in one memory block
	ku128_v unpaired;
} pcov_t;

typedef struct {
	int skip, min_pcv, min_dist, max_dist;
} pcovopt_t;

static pcov_t paircov(const rld_t *e, int len, const uint8_t *q, int skip, int max_dist, const uint64_t *sorted, hash64_t *h, uint64_t rec[3])
{
	const uint64_t mask = (uint64_t)FM_MASK30<<32 | FM_MASK30;
	fmsmem_i *iter;
	pcov_t r;
	khint_t kk;

	memset(&r, 0, sizeof(pcov_t));
	r.cov = calloc((len + 1) * 2, 1);
	r.pcv = r.cov + len + 1;
	r.len = len;
	iter = fm6_miter_init(e, len, q);
	while (fm6_miter_next(iter) >= 0) {
		int i, ret, tmp, j;
		uint64_t k, l;
		for (i = 0; i < iter->match.n; ++i) {
			fmintv_t *p = &iter->match.a[i];
			if (p->info>>63 && p->x[1] < e->mcnt[1]) { // full-length match
				tmp = p->info&FM_MASK30;
				for (j = p->info>>32&FM_MASK30; j < tmp; ++j) // update coverage
					if (r.cov[j] < 255) ++r.cov[j];
				++r.n_supp;
				if (skip <= 0 || sorted == 0) continue; // the reads are unpaired
				for (l = 0; l < p->x[2]; ++l) {
					k = sorted[p->x[1] + l] >> 2; // NB: ->x[1] corresponds to the interval of the reverse
					if ((k&1) == 0) { // reverse stand; check
						int beg, end, to_add = 0;
						kk = kh_get(64, h, k);
						if (kk != kh_end(h)) { // mate found on the forward strand
							beg = kh_val(h, kk)>>32;
							end = p->info&FM_MASK30;
							if (end - beg < max_dist) { // a proper pair
								++rec[0];
								rec[1] += end - beg;
								rec[2] += (end - beg) * (end - beg);
							} else to_add = 1; // excessive insert size
						} else to_add = 1;
						if (to_add == 1) {
							ku128_t *q;
							kv_pushp(ku128_t, r.unpaired, &q);
							q->x = k^1, q->y = p->info&mask;
							continue;
						}
						beg += skip; end -= skip;
						if (beg > end) tmp = beg, beg = end, end = tmp;
						if (beg < 0) beg = 0;
						if (end > len) end = len;
						for (j = beg; j < end; ++j)
							if (r.pcv[j] < 255) ++r.pcv[j]; // update paired coverage
						kh_del(64, h, kk);
					} else { // forward strand; add
						kk = kh_put(64, h, k^3, &ret);
						kh_val(h, kk) = p->info & mask;
					}
				}
			}
		}
	}

	for (kk = 0; kk != kh_end(h); ++kk)
		if (kh_exist(h, kk)) {
			ku128_t *q;
			kv_pushp(ku128_t, r.unpaired, &q);
			q->x = kh_key(h, kk)^2, q->y = kh_val(h, kk);
		}
	fm6_miter_destroy(iter);
	kh_clear(64, h);
	return r;
}

static void mask_pcv(int l, char *seq, const uint8_t *pcv, int skip, int min_pcv)
{
	int i, beg, end;
	for (i = 0; i < l; ++i)
		if (pcv[i] >= min_pcv) break;
	beg = i;
	if (beg == l) {
		for (i = 0; i < l; ++i)
			seq[i] = "$ACGTN"[(int)seq[i]];
		return;
	}
	for (i = 0; i < beg; ++i)
		seq[i] = beg < skip<<1? "$ACGTN"[(int)seq[i]] : "$acgtn"[(int)seq[i]];
	for (i = l - 1; i >= 0; --i)
		if (pcv[i] >= min_pcv) break;
	end = i + 1;
	for (i = end; i < l; ++i)
		seq[i] = l - end < skip<<1? "$ACGTN"[(int)seq[i]] : "$acgtn"[(int)seq[i]];
	for (i = beg; i < end; ++i)
		seq[i] = pcv[i] >= min_pcv? "$ACGTN"[(int)seq[i]] : "$acgtn"[(int)seq[i]];
}

// if unpaired, skip<=0 or sorted==0
static void paircov_all(const rld_t *e, const uint64_t *sorted, int skip, int max_dist, int n, int *len, char **s, int start, int step, int min_pcv, char *const* name,
						char *const* comment, uint64_t rec[3])
{
	int i, j;
	hash64_t *h;
	kstring_t out;

	h = kh_init(64);
	out.l = out.m = 0; out.s = 0;
	if (sorted == 0) skip = -1, min_pcv = 0; // if no rank->index map, we do not break
	for (i = start; i < n; i += step) {
		uint8_t *si = (uint8_t*)s[i];
		int l = len[i], beg, k;
		pcov_t r;
		for (j = 0; j < l; ++j)
			si[j] = seq_nt6_table[si[j]];
		if (kh_n_buckets(h) >= 256) {
			kh_destroy(64, h);
			h = kh_init(64);
		}
		r = paircov(e, l, si, skip, max_dist, sorted, h, rec);
		for (j = 0; j < l; ++j)
			r.cov[j] = r.cov[j] + 33 < 126? r.cov[j] + 33 : 126;
		if (min_pcv > 0) { // we want to break the sequence
			mask_pcv(l, (char*)si, r.pcv, skip, min_pcv);
			for (j = 0; j < l; ++j) // skip the leading lowercase letters
				if (isupper(si[j])) break;
			beg = j;
			for (j = beg + 1, k = 0; j <= l; ++j) {
				if ((islower(si[j]) || j == l) && isupper(si[j-1])) {
					out.l = 0;
					kputc('@', &out); kputs(name[i], &out); kputc('_', &out); kputw(k, &out);
					kputc('\t', &out); kputw(j - beg, &out); kputc('\t', &out); kputw(r.n_supp, &out);
					kputc('\n', &out);
					kputsn((char*)si + beg, j - beg, &out); kputsn("\n+\n", 3, &out);
					kputsn((char*)r.cov+ beg, j - beg, &out); kputc('\n', &out);
					fwrite(out.s, 1, out.l, stdout);
					++k;
				}
				if (isupper(si[j]) && islower(si[j-1])) beg = j;
			}
		} else {
			out.l = 0;
			kputc('@', &out); kputs(name[i], &out);
			if (comment[i]) {
				char *q;
				strtol(comment[i], &q, 10);
				if (q != comment[i] && isspace(*q)) {
					kputc('\t', &out); kputw(r.n_supp, &out);
					kputc('\t', &out); kputs(q + 1, &out);
				}
			}
			if (r.unpaired.n) {
				kputsn("\tUR:Z:", 6, &out);
				for (j = 0; j < r.unpaired.n; ++j) {
					kputl(r.unpaired.a[j].x, &out); kputc(',', &out);
					kputl(r.unpaired.a[j].y>>32, &out); kputc(',', &out);
					kputl(r.unpaired.a[j].y<<32>>32, &out);
					kputc(';', &out);
				}
			}
			kputc('\n', &out);
			for (j = 0; j < r.len; ++j) si[j] = "$ACGTN"[si[j]];
			kputsn((char*)si, r.len, &out); kputsn("\n+\n", 3, &out);
			kputsn((char*)r.cov, r.len, &out); kputc('\n', &out);
			fwrite(out.s, 1, out.l, stdout);
		}
		free(r.cov); free(r.unpaired.a);
	}
	kh_destroy(64, h);
	free(out.s);
}

typedef struct {
	int n, m, *l;
	char **s, **name, **comment;
} seqbuf_t;

static int fill_seqbuf(kseq_t *kseq, seqbuf_t *buf, int64_t max_len)
{
	int64_t l = 0;
	int i;
	for (i = 0; i < buf->n; ++i) {
		free(buf->s[i]);
		free(buf->name[i]);
		free(buf->comment[i]);
	}
	buf->n = 0;
	while (kseq_read(kseq) >= 0) {
		if (buf->n == buf->m) {
			buf->m = buf->m? buf->m<<1 : 256;
			buf->l = realloc(buf->l, buf->m * sizeof(int));
			buf->s = realloc(buf->s, buf->m * sizeof(void*));
			buf->name = realloc(buf->name, buf->m * sizeof(void*));
			buf->comment = realloc(buf->comment, buf->m * sizeof(void*));
		}
		buf->l[buf->n] = kseq->seq.l;
		buf->s[buf->n] = strdup(kseq->seq.s);
		buf->comment[buf->n] = kseq->comment.l? strdup(kseq->comment.s) : 0;
		buf->name[buf->n++] = strdup(kseq->name.s);
		l += kseq->seq.l;
		if (l >= max_len) break;
	}
	return buf->n;
}

typedef struct {
	const rld_t *e;
	const uint64_t *sorted;
	int start, step, skip, min_pcv, max_dist;
	seqbuf_t *buf;
	uint64_t rec[3];
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	paircov_all(w->e, w->sorted, w->skip, w->max_dist, w->buf->n, w->buf->l, w->buf->s, w->start, w->step, w->min_pcv, w->buf->name, w->buf->comment, w->rec);
	return 0;
}

int fm6_remap(const char *fn, const rld_t *e, uint64_t *sorted, int skip, int min_pcv, int max_dist, int n_threads)
{
	int i;
	kseq_t *seq;
	gzFile fp;
	worker_t *w;
	seqbuf_t *buf;
	pthread_t *tid;
	pthread_attr_t attr;
	uint64_t rec[3];
	double avg, std;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	buf = calloc(1, sizeof(seqbuf_t));
	for (i = 0; i < n_threads; ++i) {
		w[i].e = e, w[i].sorted = sorted, w[i].step = n_threads, w[i].start = i, w[i].buf = buf;
		w[i].skip = skip, w[i].min_pcv = min_pcv, w[i].max_dist = max_dist;
	}
	
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);

	while (fill_seqbuf(seq, buf, 1<<28) > 0) {
		for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], &attr, worker, w + i);
		for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	}
	rec[0] = rec[1] = rec[2] = 0;
	for (i = 0; i < n_threads; ++i)
		rec[0] += w[i].rec[0], rec[1] += w[i].rec[1], rec[2] += w[i].rec[2];
	avg = (double)rec[1] / rec[0];
	std = sqrt((double)rec[2] / rec[0] - avg * avg);
	fprintf(stderr, "[M::%s] avg = %.2f std = %.2f cap = %d\n", __func__, avg, std, (int)(avg + std * 2. + 1.499));

	free(buf->l); free(buf->s); free(buf->name); free(buf->comment); free(buf);
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
