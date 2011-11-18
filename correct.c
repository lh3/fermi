#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include "utils.h"
#include "fermi.h"
#include "rld.h"
#include "kvec.h"
#include "kstring.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint32_t)

#define SUF_LEN   8
#define SUF_SHIFT (SUF_LEN<<1)
#define SUF_NUM   (1<<SUF_SHIFT)
#define SUF_MASK  (SUF_NUM-1)

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define solid_hash(a) ((a)>>2)
#define solid_eq(a, b) ((a)>>2 == (b)>>2)
#include "khash.h"
KHASH_INIT(solid, uint32_t, uint8_t, 1, solid_hash, solid_eq)
typedef khash_t(solid) shash_t;

static volatile int g_stdout_lock;
static double g_tc, g_tr;

void seq_reverse(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

static inline double genpar_aux(double x, int64_t k)
{ // compute 1-(1-x)^k, where x<<1 and k is not so large
	int64_t i;
	double sum = 0., p = x, y = k;
	for (i = 1; i < k; ++i) {
		sum += p * y;
		p *= -x; y *= (double)(k - i) / (i + 1);
		if (p * y / sum < 1e-6) break;
	}
	return sum;
}

void fm_ec_genpar(int64_t n, int l, double cov, double p, int *_w, int *_T)
{
	int w, k;
	int64_t L;
	double e, qc, qe;
	L = (int64_t)((double)n * l / cov + .499);
	e = genpar_aux(p, l) * n;
	for (w = 8; w < l; ++w) {
		double q, D;
		q = genpar_aux(p, w) * (1 - p) * genpar_aux(pow(.25, w), L) * .75;
		D = genpar_aux(q, l - w) * pow(1 - p, l) * n;
		if (D < 0.0001 * e) break;
	}
	qc = (double)(l - w) / L * pow(1 - p, w + 1);
	qe = (double)(l - w) / L * (1./3.) * p * pow(1 - p, w);
	for (k = 1; k < (int)cov + 1; ++k)
		if (pow(qc, k) * pow(1 - qc, n - k) > pow(qe, k) * pow(1 - qe, n - k)) break;
	k += 2;
	fprintf(stderr, "[M::%s] HiTEC parameters for n=%ld, l=%d and c=%.1f: w_M=%d, T(w_M)=%d\n", __func__, (long)n, l, cov, w, k);
	*_w = w; *_T = k;
}

static void ec_collect(const rld_t *e, const fmecopt_t *opt, int len, const uint8_t *seq, shash_t *solid, int64_t cnt[2])
{
	int i, ret, shift = (opt->w - len - 1) * 2;
	kstring_t str;
	fmintv_v stack;
	fmintv_t ok[6], ik;

	assert(len > 0 && opt->w > len);
	kv_init(stack);
	str.m = opt->w + 1;
	str.s = calloc(str.m, 1);
	str.l = 0;
	fm6_set_intv(e, seq[0], ik); // descend to the root of the subtree to be processed
	for (i = 1; i < len; ++i) {
		fm6_extend(e, &ik, ok, 1);
		ik = ok[(int)seq[i]];
	}

	ik.info = len<<4;
	kv_push(fmintv_t, stack, ik);
	while (stack.n) {
		int c;
		ik = kv_pop(stack);
		fm6_extend(e, &ik, ok, 1);
		str.l = (ik.info>>4) - len;
		if (str.l) str.s[str.l - 1] = ik.info&0xf;
		if (ik.info>>4 == opt->w) { // keep the k-mer
			uint32_t key;
			int max_c;
			khint_t k;
			uint64_t max, rest;
			double r;
			for (c = 1, max = 0, max_c = 6; c <= 4; ++c)
				if (ok[c].x[2] > max)
					max = ok[c].x[2], max_c = c;
			if (max < opt->T) continue;
			++cnt[0];
			rest = ik.x[2] - max;
			r = rest == 0? 31. : (double)max / rest;
			if (r > 31.) r = 31.;
			if (rest <= 7 && r >= opt->T) ++cnt[1];
			for (i = 0, key = 0; i < str.l; ++i) key = (uint32_t)str.s[i]<<shift | key>>2;
			key = key<<2 | (max_c - 1);
			k = kh_put(solid, solid, key, &ret);
			kh_val(solid, k) = rest > 7? 0 : ((int)(r + .499)) << 3 | rest; // if more than 7, do not make corrections in future
			//kh_val(solid, k) = (int)(r + .499) << 3 | (rest < 7? rest : 7);
		} else { // descend
			for (c = 4; c >= 1; --c) { // FIXME: ambiguous bases are skipped
				if (ok[c].x[2] >= opt->T + 1) {
					ok[c].info = ((ik.info>>4) + 1)<<4 | (c - 1);
					kv_push(fmintv_t, stack, ok[c]);
				}
			}
		}
	}

	free(stack.a); free(str.s);
}

static int ec_fix1(const fmecopt_t *opt, shash_t *const* solid, kstring_t *s)
{
	int i, l, cnt = 0, shift = (opt->w - 1) * 2, mod = 0;
	uint64_t x;
	if (s->l <= opt->w) return -1; // correction failed
	for (i = s->l - 1, l = 0, x = 0; i >= 0; --i) {
		if (s->s[i] == 5) {
			if (mod == 0) {
				x = 0, l = 0;
				continue;
			} else s->s[i] = mod;
		}
		x = (uint64_t)(s->s[i]-1)<<shift | x>>2;
		if (++l >= opt->w) {
			khint_t k;
			const shash_t *h = solid[x & SUF_MASK];
			k = kh_get(solid, h, x>>SUF_SHIFT<<2);
//			if (fm_verbose >= 10) fprintf(stderr, "%d\t%c\t%d\t%c\t%c\n", i, "$ACGTN"[s->s[i]], k == kh_end(h)?-1:kh_val(h,k)>>3, "$ACGTN"[mod], k==kh_end(h)?'?':"ACGTN"[kh_key(h, k)&3]);
			if (k == kh_end(h) && mod) {
				s->s[i] = mod;
				x = (x & ((1LLU<<shift) - 1)) | (uint64_t)(mod - 1) << shift;
				k = kh_get(solid, h, x>>SUF_SHIFT<<2);
				++cnt;
			}
			mod = (k != kh_end(h) && kh_val(h, k)>>3 >= opt->min_ratio && i && s->s[i-1] != (kh_key(h, k)&3) + 1)? (kh_key(h, k)&3) + 1 : 0;
		} else mod = 0;
	}
	return cnt;
}

static void ec_fix(const rld_t *e, const fmecopt_t *opt, shash_t *const* solid, int n_seqs, char **seq, char **qual)
{
	extern unsigned char seq_nt6_table[];
	int i, j, cnt[2];
	kstring_t str;

	str.s = 0; str.l = str.m = 0;
	for (i = 0; i < n_seqs; ++i) {
		str.l = 0;
		kputs(seq[i], &str);
		for (j = 0; j < str.l; ++j)
			str.s[j] = seq_nt6_table[(int)str.s[j]];
		seq_revcomp6(str.l, (uint8_t*)str.s);
		cnt[0] = ec_fix1(opt, solid, &str);
		if (cnt[0] < 0) continue;
		seq_revcomp6(str.l, (uint8_t*)str.s);
		cnt[1] = ec_fix1(opt, solid, &str);
		for (j = 0; j < str.l; ++j)
			seq[i][j] = seq_nt6_table[(int)seq[i][j]] == str.s[j]? toupper(seq[i][j]) : "$acgtn"[(int)str.s[j]];
	}
	free(str.s);
}

typedef struct {
	const rld_t *e;
	const fmecopt_t *opt;
	shash_t **solid;
	int64_t cnt[2];
	int n_seqs, tid;
	uint32_t *seqs;
} worker1_t;

static void *worker1(void *data)
{
	worker1_t *w = (worker1_t*)data;
	int i, j;
	for (i = 0; i < w->n_seqs; ++i) {
		uint8_t seq[SUF_LEN+1];
		for (j = 0; j < SUF_LEN; ++j)
			seq[j] = (w->seqs[i]>>(j*2)&0x3) + 1;
		ec_collect(w->e, w->opt, SUF_LEN, seq, w->solid[i], w->cnt);
	}
	return 0;
}

#define BATCH_SIZE 1000000

typedef struct {
	const rld_t *e;
	const fmecopt_t *opt;
	shash_t *const* solid;
	int n_seqs;
	char **seq, **qual;
	uint64_t *id; // this wastes memory, but should not be a big deal
} worker2_t;

static void *worker2(void *data)
{
	worker2_t *w = (worker2_t*)data;
	ec_fix(w->e, w->opt, w->solid, w->n_seqs, w->seq, w->qual);
	return 0;
}

int fm6_ec_correct(const rld_t *e, const fmecopt_t *opt, const char *fn, int n_threads)
{
	int j;
	int64_t i, cnt[2];
	shash_t **solid;
	pthread_t *tid;
	pthread_attr_t attr;

	if (opt->w <= SUF_LEN) {
		fprintf(stderr, "[E::%s] excessively small `-w'. Please increase manually to at least %d.\n", __func__, SUF_LEN + 1);
		return 1;
	}
	// initialize "solid" and "tid"
	assert(n_threads <= SUF_NUM);
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	solid = calloc(SUF_NUM, sizeof(void*));
	for (j = 0; j < SUF_NUM; ++j) solid[j] = kh_init(solid);
	cnt[0] = cnt[1] = 0;

	{ // initialize and launch worker1
		worker1_t *w1;
		int max_seqs;
		g_tc = cputime(); g_tr = realtime();
		w1 = calloc(n_threads, sizeof(worker1_t));
		max_seqs = (SUF_NUM + n_threads - 1) / n_threads;
		for (j = 0; j < n_threads; ++j) {
			w1[j].seqs = calloc(max_seqs, 4);
			w1[j].solid = calloc(max_seqs, sizeof(void*));
			w1[j].e = e, w1[j].opt = opt, w1[j].tid = j;
		}
		for (i = 0, j = 0; i < SUF_NUM; ++i) { // assign seqs and hash tables
			w1[j].solid[w1[j].n_seqs] = solid[i];
			w1[j].seqs[w1[j].n_seqs++] = i;
			if (++j == n_threads) j = 0;
		}
		for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker1, w1 + j);
		for (j = 0; j < n_threads; ++j) {
			pthread_join(tid[j], 0);
			free(w1[j].seqs); free(w1[j].solid);
			cnt[0] += w1[j].cnt[0], cnt[1] += w1[j].cnt[1];
		}
		free(w1);
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] collected %ld informative and %ld ambiguous k-mers in %.3f sec (%.3f wall clock)\n",
					__func__, (long)cnt[1], (long)(cnt[0] - cnt[1]), cputime() - g_tc, realtime() - g_tr);
	}

	{ // initialize and launch worker2
		gzFile fp;
		kseq_t *seq;
		worker2_t *w2;
		int max_seqs;
		uint64_t k, id = 0, pre_id = 0;
		kstring_t out;

		g_tc = cputime(); g_tr = realtime();
		out.m = out.l = 0; out.s = 0;
		fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
		seq = kseq_init(fp);
		w2 = calloc(n_threads, sizeof(worker2_t));
		max_seqs = ((BATCH_SIZE < e->mcnt[1]/2? BATCH_SIZE : e->mcnt[1]/2) + n_threads - 1) / n_threads;
		for (j = 0; j < n_threads; ++j) {
			w2[j].e = e, w2[j].solid = solid, w2[j].opt = opt;
			w2[j].seq  = calloc(max_seqs, sizeof(void*));
			w2[j].qual = calloc(max_seqs, sizeof(void*));
			w2[j].id   = calloc(max_seqs, 8);
		}
		for (;;) {
			int ret;
			worker2_t *w;
			ret = kseq_read(seq);
			if (ret < 0 || (id && id%BATCH_SIZE == 0)) {
				for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker2, w2 + j);
				for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
				for (j = 0; j < n_threads; ++j) w2[j].n_seqs = 0;
				for (k = pre_id; k < id; ++k) {
					w = &w2[k%n_threads];
					out.l = 0;
					kputc('@', &out); kputw(w->id[w->n_seqs]>>1, &out); kputc('\n', &out);
					kputs(w->seq[w->n_seqs], &out);
					kputsn("\n+\n", 3, &out); kputs(w->qual[w->n_seqs], &out);
					puts(out.s);
					free(w->seq[w->n_seqs]); free(w->qual[w->n_seqs]);
					++w->n_seqs;
				}
				for (j = 0; j < n_threads; ++j) w2[j].n_seqs = 0;
				if (fm_verbose >= 3)
					fprintf(stderr, "[M::%s] corrected errors in %ld reads in %.3f CPU seconds (%.3f wall clock)\n",
							__func__, (long)id, cputime() - g_tc, realtime() - g_tr);
				pre_id = id;
			}
			if (ret < 0) break;
			w = &w2[id%n_threads];
			w->seq[w->n_seqs] = strdup(seq->seq.s);
			if (seq->qual.l == 0) { // if no quality, set to 20
				w->qual[w->n_seqs] = malloc(seq->seq.l + 1);
				for (j = 0; j < seq->seq.l; ++j)
					w->qual[w->n_seqs][j] = 33 + 15;
				w->qual[w->n_seqs][j] = 0;
			} else w->qual[w->n_seqs] = strdup(seq->qual.s);
			w->id[w->n_seqs++] = id<<1;
			++id;
		}
		for (j = 0; j < n_threads; ++j) {
			free(w2[j].seq); free(w2[j].qual); free(w2[j].id);
		}
		free(w2); free(out.s);
		kseq_destroy(seq);
		gzclose(fp);
	}

	// free
	for (j = 0; j < SUF_NUM; ++j) kh_destroy(solid, solid[j]);
	free(solid); free(tid);
	return 0;
}
