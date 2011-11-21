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

#define SUF_LEN   8
#define SUF_SHIFT (SUF_LEN<<1)
#define SUF_NUM   (1<<SUF_SHIFT)
#define SUF_MASK  (SUF_NUM-1)
#define MAX_KMER  (SUF_LEN + 15)

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
void ks_introsort_128y(size_t n, fm128_t *a); // in msg.c
void ks_heapup_128y(size_t n, fm128_t *a);
void ks_heapdown_128y(size_t i, size_t n, fm128_t *a);

/****************************
 * Compute HiTEC parameters *
 ****************************/

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
	if (fm_verbose >= 3) fprintf(stderr, "[M::%s] HiTEC parameters for n=%ld, l=%d and c=%.1f: w_M=%d, T(w_M)=%d\n", __func__, (long)n, l, cov, w, k);
	if (k > MAX_KMER) {
		k = MAX_KMER;
		if (fm_verbose >= 2) fprintf(stderr, "[W::%s] Set k-mer length to the maximum length %d\n", __func__, MAX_KMER);
	}
	*_w = w; *_T = k;
}

/***********************
 * Collect good k-mers *
 ***********************/

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
			if (max < opt->min_occ) continue; // then in the following max_c<6
			++cnt[0];
			rest = ik.x[2] - max - ok[0].x[2] - ok[5].x[2];
			r = rest == 0? max : (double)max / rest;
			if (r > 31.) r = 31.; // we have maximally 5 bits of information (i.e. [0,31])
			if (rest <= 7 && r >= opt->min_occ) ++cnt[1];
			for (i = 0, key = 0; i < str.l; ++i)
				key = (uint32_t)str.s[i]<<shift | key>>2;
			key = key<<2 | (max_c - 1);
			k = kh_put(solid, solid, key, &ret);
			kh_val(solid, k) = (int)(r + .499) << 3 | (rest < 7? rest : 7);
		} else { // descend
			for (c = 4; c >= 1; --c) { // ambiguous bases are skipped
				if (ok[c].x[2] >= opt->min_occ) {
					ok[c].info = ((ik.info>>4) + 1)<<4 | (c - 1);
					kv_push(fmintv_t, stack, ok[c]);
				}
			}
		}
	}

	free(stack.a); free(str.s);
}

/******************
 * Correct errors *
 ******************/

typedef struct {
	fm128_v heap;
	fm32_v stack;
} fixaux_t;

static inline void save_state(fixaux_t *fa, const fm128_t *p, int c, int score, int shift)
{
	fm128_t w;
	if (score < 0) score = 0;
	w.x = (uint64_t)c<<shift | p->x>>2;
	w.y = (uint64_t)((p->y>>48) + score)<<48 | fa->stack.n<<16 | ((p->y&0xffff) - 1);
	kv_push(uint32_t, fa->stack, c<<28 | (uint32_t)(p->y>>16));
	kv_push(fm128_t, fa->heap, w);
	ks_heapup_128y(fa->heap.n, fa->heap.a);
}

#define TIMES_FACTOR 10
#define DIFF_FACTOR  13
#define MISS_PENALTY 40
#define MAX_HEAP     64

static int ec_fix1(const fmecopt_t *opt, shash_t *const* solid, kstring_t *s, char *qual, fixaux_t *fa)
{
	int i, l, shift = (opt->w - 1) << 1, n_rst = 0, first = 0, qsum;
	fm128_t z, rst[2];

	if (s->l <= opt->w) return -1;
	// get the initial k-mer
	fa->heap.n = fa->stack.n = 0;
	z.x = z.y = 0;
	for (i = s->l - 1, l = 0; i > 0 && l < opt->w; --i)
		if (s->s[i] == 5) z.x = 0, l = 0;
		else z.x = (uint64_t)(s->s[i]-1)<<shift | z.x>>2, ++l;
	if (i == 0) return 0xffff; // no good k-mer
	// the first element in the heap
	kv_push(uint32_t, fa->stack, 0);
	z.y = i + 1;
	kv_push(fm128_t, fa->heap, z);
	// traverse
	while (fa->heap.n) {
		const shash_t *h;
		khint_t k;
		int parent;
		// get the best so far
		z = fa->heap.a[0];
		fa->heap.a[0] = kv_pop(fa->heap);
		ks_heapdown_128y(0, fa->heap.n, fa->heap.a);
		if ((z.y&0xffff) == 0) {
			rst[n_rst++] = z;
			if (n_rst == 2) break;
			continue;
		}
		i = (z.y&0xffff) - 1; parent = (int)(z.y>>16);
		// check the hash table
		h = solid[z.x & SUF_MASK];
		k = kh_get(solid, h, z.x>>SUF_SHIFT<<2);
		if (k != kh_end(h)) { // this (k+1)-mer has more than opt->min_occ occurrences
			first = 1;
			if (s->s[i] != (kh_key(h, k)&3) + 1) { // the read base is different from the best base
				int tmp, score, max = (kh_val(h, k)&7)? (kh_val(h, k)&7) * (kh_val(h, k)>>3) : kh_val(h, k)>>3;
				score = (max - (kh_val(h, k)&7)) * DIFF_FACTOR; // score for the best stack path
				tmp = (kh_val(h, k)&7)? (kh_val(h, k)>>3) * TIMES_FACTOR : 10000;
				if (tmp < score) score = tmp;
				tmp = (7 - (kh_val(h, k)&7)) * DIFF_FACTOR;
				if (tmp < score) score = tmp;
				// if we have too many possibilities, keep the better path among the two
				if (s->s[i] != 5 && (fa->heap.n + 2 <= MAX_HEAP || score < qual[i]-33))
					save_state(fa, &z, s->s[i] - 1, score, shift); // the read path
				if (s->s[i] == 5 || fa->heap.n + 2 <= MAX_HEAP || score > qual[i]-33)
					save_state(fa, &z, kh_key(h, k)&3, qual[i]-33, shift); // the stack path
			} else save_state(fa, &z, s->s[i] - 1, 0, shift);
		} else {
			int score = first? MISS_PENALTY - (qual[i] - 33) : 0;
			save_state(fa, &z, s->s[i] - 1, score > 0? score : 0, shift);
		}
	}
	assert(n_rst == 1 || n_rst == 2);
	if (rst[0].y>>48 == 0) return 0x10000000; // no corrections are made
	// backtrack
	i = 0; qsum = 0; l = (uint32_t)(rst[0].y>>16);
	while (l) {
		if (s->s[i] - 1 != fa->stack.a[l]>>28) {
			s->s[i] = (fa->stack.a[l]>>28) + 1;
			qsum += qual[i] - 33;
			qual[i] = 34; // reduce the base quality
		}
		++i;
		l = fa->stack.a[l]<<4>>4;
	}
	l = n_rst == 2? (rst[1].y>>48) - (rst[0].y>>48) : 0xfff;
	if (l > 0xfff) l = 0xfff;
	return qsum | l<<18;
}

static void ec_fix(const rld_t *e, const fmecopt_t *opt, shash_t *const* solid, int n_seqs, char **seq, char **qual, int *info)
{
	extern unsigned char seq_nt6_table[];
	int i, j, ret0, ret1;
	kstring_t str;
	fixaux_t fa;

	memset(&fa, 0, sizeof(fixaux_t));
	str.s = 0; str.l = str.m = 0;
	for (i = 0; i < n_seqs; ++i) {
		str.l = 0;
		kputs(seq[i], &str);
		for (j = 0; j < str.l; ++j)
			str.s[j] = seq_nt6_table[(int)str.s[j]];
		seq_revcomp6(str.l, (uint8_t*)str.s); // to the reverse complement strand
		seq_reverse(str.l, (uint8_t*)qual[i]);
		ret0 = ec_fix1(opt, solid, &str, qual[i], &fa); // 0x7fff0000 if no correction; 0xffff if too short
		seq_reverse(str.l, (uint8_t*)qual[i]); // back to the forward strand
		seq_revcomp6(str.l, (uint8_t*)str.s);
		if (ret0 != 0xffff) { // then we need to correct in the reverse direction
			uint64_t k, l;
			ret1 = ec_fix1(opt, solid, &str, qual[i], &fa);
			info[i] = ((ret0&0xffff) + (ret1&0xffff)) | (ret0>>18 < ret1>>18? ret0>>18 : ret1>>18)<<18;
			// FIXME: the following can be accelerated because we may know if the last k-mer has a match
			for (j = 0; j < opt->w; ++j)
				if (str.s[j] == 5) break;
			if (j == opt->w) {
				fm_backward_search(e, opt->w, (uint8_t*)str.s, &k, &l);
				if (l - k > 1) info[i] |= 1<<16;
			}
			for (j = str.l - opt->w; j < str.l; ++j)
				if (str.s[j] == 5) break;
			if (j == str.l) {
				fm_backward_search(e, opt->w, (uint8_t*)str.s + (str.l - opt->w), &k, &l);
				if (l - k > 1) info[i] |= 1<<17;
			}
		} else info[i] = ret0;
		for (j = 0; j < str.l; ++j)
			seq[i][j] = seq_nt6_table[(int)seq[i][j]] == str.s[j]? toupper(seq[i][j]) : "$acgtn"[(int)str.s[j]];
	}
	free(str.s);
	free(fa.heap.a); free(fa.stack.a);
}

/************************
 * Multi-thread workers *
 ************************/

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
	int n_seqs, *info;
	char **seq, **qual;
	uint64_t *id; // this wastes memory, but should not be a big deal
} worker2_t;

static void *worker2(void *data)
{
	worker2_t *w = (worker2_t*)data;
	ec_fix(w->e, w->opt, w->solid, w->n_seqs, w->seq, w->qual, w->info);
	return 0;
}

/******************
 * The key portal *
 ******************/

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
			w2[j].info = calloc(max_seqs, sizeof(int));
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
					kputc('@', &out); kputl(w->id[w->n_seqs]>>1, &out);
					kputc('_', &out); kputw(w->info[w->n_seqs]>>16&3, &out);
					kputc('_', &out); kputw(w->info[w->n_seqs]&0xffff, &out);
					kputc('_', &out); kputw(w->info[w->n_seqs]>>18, &out); kputc('\n', &out);
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
			free(w2[j].seq); free(w2[j].qual); free(w2[j].id); free(w2[j].info);
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
