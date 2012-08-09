#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include "priv.h"
#include "kvec.h"
#include "kstring.h"

static int SUF_LEN, SUF_NUM;

#include <zlib.h>
#include "kseq.h"
KSEQ_DECLARE(gzFile)

#define solid_hash(a) ((a)>>2)
#define solid_eq(a, b) ((a)>>2 == (b)>>2)
#include "khash.h"
KHASH_INIT(solid, uint32_t, uint8_t, 1, solid_hash, solid_eq)
typedef khash_t(solid) shash_t;

static double g_tc, g_tr;

static void compute_SUF(int suf_len)
{
	SUF_LEN   = suf_len;
	SUF_NUM   = 1<<(SUF_LEN<<1);
}

/***********************
 * Collect good k-mers *
 ***********************/

static void ec_collect(const rld_t *e, const fmecopt_t *opt, int len, const fmintv_t *suf_intv, shash_t *solid, int64_t cnt[2])
{
	int i, ret, shift = (opt->w - len - 1) * 2;
	kstring_t str;
	fmintv_v stack;
	fmintv_t ok[6], ik;

	if (suf_intv->x[2] == 0) return;
	assert(len > 0 && opt->w > len);
	kv_init(stack);
	str.m = opt->w + 1; str.l = 0;
	str.s = calloc(str.m, 1);
	ik = *suf_intv;
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
	ku128_v heap;
	fm64_v stack;
} fixaux_t;

static inline void save_state(fixaux_t *fa, const ku128_t *p, int c, int score, int shift, int has_match)
{
	ku128_t w;
	if (score < 0) score = 0;
	if (c >= 4) c = 0;
	w.x = (uint64_t)c<<shift | p->x>>2;
	// the structure of w.y - score:16, pos_in_stack:32, seq_pos:16
	w.y = (uint64_t)((p->y>>48) + score)<<48 | fa->stack.n<<16 | ((p->y&0xffff) - 1);
	// structure of a stack element - read_pos:32, base:3, has_match:1, parent_pos_in_stack:28
	kv_push(uint64_t, fa->stack, ((p->y&0xffff) - 1)<<32 | (uint32_t)c<<29 | has_match<<28 | (uint32_t)(p->y>>16));
	kv_push(ku128_t, fa->heap, w);
	ks_heapup_128y(fa->heap.n, fa->heap.a);
}

#define RATIO_FACTOR 10
#define DIFF_FACTOR  13
#define MAX_HEAP    256
#define MAX_SC_DIFF  60
#define MAX_QUAL     40
#define MISS_PENALTY 10
#define MIN_OCC       5
#define MIN_OCC_RATIO 0.8

static int ec_fix1(const fmecopt_t *opt, shash_t *const* solid, kstring_t *s, char *qual, fixaux_t *fa, uint64_t *n_query)
{
	int i, q, l, shift = (opt->w - 1) << 1, n_rst = 0, qsum, no_hits = 1, score_diff;
	ku128_t z, rst[2];

	if (s->l <= opt->w) return 0xffff;
	// get the initial k-mer
	fa->heap.n = fa->stack.n = 0;
	z.x = z.y = 0;
	for (i = s->l - 1, l = 0; i > 0 && l < opt->w; --i)
		if (s->s[i] == 5) z.x = 0, l = 0;
		else z.x = (uint64_t)(s->s[i]-1)<<shift | z.x>>2, ++l;
	if (i == 0) return 0xffff; // no good k-mer
	// the first element in the heap
	kv_push(uint64_t, fa->stack, 0);
	z.y = i + 1;
	kv_push(ku128_t, fa->heap, z);
	// traverse
	while (fa->heap.n) {
		const shash_t *h;
		khint_t k;
		// get the best so far
		z = fa->heap.a[0];
		fa->heap.a[0] = kv_pop(fa->heap);
		ks_heapdown_128y(0, fa->heap.n, fa->heap.a);
		if ((z.y&0xffff) == 0) {
			rst[n_rst++] = z;
			if (n_rst == 2) break;
			continue;
		}
		if (n_rst && (int)(z.y>>48) > (int)(rst[0].y>>48) + MAX_SC_DIFF) break;
		i = (z.y&0xffff) - 1;
		q = qual[i] - 33 < MAX_QUAL? qual[i] - 33 : MAX_QUAL;
		if (q < 3) q = 3;
		// check the hash table
		h = solid[z.x & (SUF_NUM - 1)];
		k = kh_get(solid, h, z.x>>(SUF_LEN<<1)<<2);
		++*n_query;
		if (k != kh_end(h)) { // this (k+1)-mer has more than opt->min_occ occurrences
			no_hits = 0;
			if (s->s[i] != (kh_key(h, k)&3) + 1) { // the read base is different from the best base
				int v = kh_val(h, k); // recall that v is packed as - "(best_depth/rest_depth)<<3 | rest_depth" or "best_detph<<3 | 0"
				int tmp, penalty, max = (v&7)? (v&7) * (v>>3) : v>>3; // max is the approximate depth of the best base
				// compute the penalty for the best stack path
				penalty = (max - (v&7)) * DIFF_FACTOR;
				if (max - (v&7) < 1) penalty = 1;
				tmp = (v&7)? (v>>3) * RATIO_FACTOR : 10000;
				if (tmp < penalty) penalty = tmp;
				tmp = (7 - (v&7)) * DIFF_FACTOR;
				if (tmp < penalty) penalty = tmp;
				if (penalty < 1) penalty = 1;
				// if we have too many possibilities, keep the better path among the two
				if (s->s[i] != 5 && (fa->heap.n + 2 <= MAX_HEAP || penalty < q))
					save_state(fa, &z, s->s[i] - 1, penalty, shift, 1); // the read path
				if (s->s[i] == 5 || fa->heap.n + 2 <= MAX_HEAP || penalty > q)
					save_state(fa, &z, kh_key(h, k)&3, q, shift, 1); // the stack path
			} else { // the read base is the same as the best base
				ku128_t z0 = z;
				int i0 = i;
				int v = kh_val(h, k), occ_last = (v&7)? (v&7) * ((v>>3)+1) : v>>3;
				if ((v&7) <= 0 && opt->step > 1) {
					while (i0 > 0) {
						for (i = (z.y&0xffff) - 1, l = 0; i >= 1 && l < opt->step && s->s[i] < 5; --i, ++l)
							z.x = (uint64_t)(s->s[i]-1)<<shift | z.x>>2; // look opt->w/2 mer ahead
						if (s->s[i] == 5) break;
						h = solid[z.x & (SUF_NUM - 1)];
						k = kh_get(solid, h, z.x>>(SUF_LEN<<1)<<2);
						++*n_query;
						if (k != kh_end(h) && s->s[i] == (kh_key(h, k)&3) + 1) { // in the hash table and the read base is the best
							int v = kh_val(h, k), occ = (v&7)? (v&7) * ((v>>3)+1) : v>>3; // occ is the occurrences of the k-mer
							if ((v&7) <= 1 && occ >= MIN_OCC && (double)occ / occ_last >= MIN_OCC_RATIO) { // if occ is good enough, jump again
								z.y = z.y>>16<<16 | (i + 1);
								z0 = z; i0 = i;
								occ_last = occ;
							} else break; // if not good, reject and stop jumping
						} else break;
					}
				}
				save_state(fa, &z0, s->s[i0] - 1, 0, shift, 1);
			}
		} else save_state(fa, &z, s->s[i] - 1, MISS_PENALTY + (MAX_QUAL - q), shift, 0);
	}
	assert(n_rst == 1 || n_rst == 2);
	score_diff = n_rst == 1? MAX_SC_DIFF : (int)(rst[1].y>>48) - (int)(rst[0].y>>48);
	assert(score_diff >= 0);
	if (score_diff >= MAX_SC_DIFF) score_diff = MAX_SC_DIFF;
	if (rst[0].y>>48 == 0) return score_diff<<18; // no corrections are made
	// backtrack
	qsum = 0; l = (uint32_t)(rst[0].y>>16);
	while (l) {
		i = fa->stack.a[l]>>32;
		if (s->s[i] - 1 != (uint32_t)fa->stack.a[l]>>29) {
			s->s[i] = ((uint32_t)fa->stack.a[l]>>29) + 1;
			qsum += qual[i] - 33;
		} else if (((uint32_t)fa->stack.a[l]>>28&1) && qual[i] < 37) qual[i] = 37;
		l = (uint32_t)fa->stack.a[l]<<4>>4;
	}
	// return value: score_diff:14, (empty):2, sum_modified_qual:16
	return qsum | score_diff<<18 | no_hits<<17;
}

static uint64_t ec_fix(const rld_t *e, const fmecopt_t *opt, shash_t *const* solid, int n_seqs, char **seq, char **qual, int *info)
{
	int i, j, ret0, ret1, n_lower;
	uint64_t n_query = 0;
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
		ret0 = ec_fix1(opt, solid, &str, qual[i], &fa, &n_query); // 0x7fff0000 if no correction; 0xffff if too short
		seq_reverse(str.l, (uint8_t*)qual[i]); // back to the forward strand
		seq_revcomp6(str.l, (uint8_t*)str.s);
		if (ret0 != 0xffff) { // then we need to correct in the reverse direction
			ret1 = ec_fix1(opt, solid, &str, qual[i], &fa, &n_query);
			info[i] = ((ret0&0xffff) + (ret1&0xffff)) | (ret0>>18 < ret1>>18? ret0>>18 : ret1>>18)<<18;
			if ((ret0>>17&1) && (ret1>>17&1)) info[i] |= 1<<16;
		} else info[i] = ret0;
		for (j = 0, n_lower = 0; j < str.l; ++j) {
			seq[i][j] = seq_nt6_table[(int)seq[i][j]] == str.s[j]? toupper(seq[i][j]) : "$acgtn"[(int)str.s[j]];
			if (islower(seq[i][j])) ++n_lower, qual[i][j] = 36;
		}
		if ((double)n_lower / str.l > opt->max_corr) info[i] |= 1<<16;
		if (info[i]>>18 <= 10) info[i] |= 1<<16;
	}
	free(str.s);
	free(fa.heap.a); free(fa.stack.a);
	return n_query;
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
	const fmintv_t *top;
} worker1_t;

static void *worker1(void *data)
{
	worker1_t *w = (worker1_t*)data;
	int i;
	for (i = 0; i < w->n_seqs; ++i)
		ec_collect(w->e, w->opt, SUF_LEN, &w->top[w->seqs[i]], w->solid[i], w->cnt);
	return 0;
}

#define BATCH_SIZE 1000000

typedef struct {
	const rld_t *e;
	const fmecopt_t *opt;
	shash_t *const* solid;
	int n_seqs, *info;
	char **seq, **qual;
	uint64_t n_query;
} worker2_t;

static void *worker2(void *data)
{
	worker2_t *w = (worker2_t*)data;
	w->n_query = ec_fix(w->e, w->opt, w->solid, w->n_seqs, w->seq, w->qual, w->info);
	return 0;
}

/******************
 * The key portal *
 ******************/

#define MAX_KMER 27

int fm6_ec_correct(const rld_t *e, fmecopt_t *opt, const char *fn, int _n_threads)
{
	int j, n_threads;
	int64_t i, cnt[2];
	shash_t **solid;
	pthread_t *tid;
	pthread_attr_t attr;

	if (opt->w < 0) { // determine k-mer
		opt->w = (int)(log(e->mcnt[0]) / log(4) + 8.499);
		if (opt->w >= MAX_KMER) opt->w = MAX_KMER;
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] set k-mer length to %d\n", __func__, opt->w);
	}
	compute_SUF(opt->w > 15? opt->w - 15 : 1);
	// initialize "solid" and "tid"
	assert(_n_threads <= SUF_NUM);
	tid = (pthread_t*)calloc(_n_threads, sizeof(pthread_t));
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	solid = calloc(SUF_NUM, sizeof(void*));
	for (j = 0; j < SUF_NUM; ++j) solid[j] = kh_init(solid);
	cnt[0] = cnt[1] = 0;

	{ // initialize and launch worker1
		worker1_t *w1;
		fmintv_t *top;
		int max_seqs;
		n_threads = _n_threads%2? _n_threads : _n_threads - 1;
		g_tc = cputime(); g_tr = realtime();
		w1 = calloc(n_threads, sizeof(worker1_t));
		max_seqs = (SUF_NUM + n_threads - 1) / n_threads;
		top = fm6_traverse(e, SUF_LEN);
		if (fm_verbose >= 3)
			fprintf(stderr, "[M::%s] traverse the trie up to depth %d in %.3f sec\n", __func__, SUF_LEN, cputime() - g_tc);
		g_tc = cputime(); g_tr = realtime();
		for (j = 0; j < n_threads; ++j) {
			w1[j].seqs = calloc(max_seqs, 4);
			w1[j].solid = calloc(max_seqs, sizeof(void*));
			w1[j].e = e, w1[j].top = top, w1[j].opt = opt, w1[j].tid = j;
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
		free(top);
	}

	{ // initialize and launch worker2
		gzFile fp;
		kseq_t *seq;
		worker2_t *w2;
		int max_seqs;
		uint64_t k, id = 0, pre_id = 0;
		kstring_t out;

		n_threads = _n_threads;
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
		}
		for (;;) {
			int ret;
			worker2_t *w;
			ret = kseq_read(seq);
			if (ret < 0 || (id && id%BATCH_SIZE == 0)) {
				uint64_t n_query = 0, n_seq = 0;
				for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker2, w2 + j);
				for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
				for (j = 0; j < n_threads; ++j) {
					n_seq += w2[j].n_seqs;
					n_query += w2[j].n_query;
					w2[j].n_seqs = w2[j].n_query = 0;
				}
				if (fm_verbose >= 3)
					fprintf(stderr, "[M::%s] corrected errors in %ld reads in %.3f CPU seconds (%.3f wall clock); %.2f lookups per read\n",
							__func__, (long)id, cputime() - g_tc, realtime() - g_tr, (double)n_query / n_seq);
				for (k = pre_id; k < id; ++k) {
					int is_bad = 0;
					w = &w2[k%n_threads];
					if (opt->is_paired) {
						if (w->info[w->n_seqs]>>16&1) is_bad = 1;
						else if (k&1) { // second in a pair
							worker2_t *m = &w2[(k^1)%n_threads];
							if (m->info[m->n_seqs - 1]>>16&1) is_bad = 1;
						} else { // first in a pair
							worker2_t *m = &w2[(k^1)%n_threads];
							if (m->info[m->n_seqs + (m == w)]>>16&1) is_bad = 1; // FIXED on 20111206: special treatment for single-thread
						}
					} else if (w->info[w->n_seqs]>>16&1) is_bad = 1;
					if (!is_bad || opt->keep_bad) {
						int tmp = 0;
						out.l = 0;
						kputc('@', &out); kputl(opt->is_paired? k>>1 : k, &out);
						kputc(opt->is_paired? ' ':'_', &out); kputw(w->info[w->n_seqs]&0xffff, &out);
						kputc(opt->is_paired? ' ':'_', &out); kputw(w->info[w->n_seqs]>>18, &out); kputc('\n', &out);
						tmp = strlen(w->seq[w->n_seqs]);
						if (opt->trim_l && opt->trim_l < tmp) tmp = opt->trim_l;
						kputsn(w->seq[w->n_seqs], tmp, &out);
						kputsn("\n+\n", 3, &out); kputsn(w->qual[w->n_seqs], tmp, &out);
						puts(out.s);
					}
					free(w->seq[w->n_seqs]); free(w->qual[w->n_seqs]);
					++w->n_seqs;
				}
				for (j = 0; j < n_threads; ++j) w2[j].n_seqs = 0;
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
			++w->n_seqs;
			++id;
		}
		for (j = 0; j < n_threads; ++j) {
			free(w2[j].seq); free(w2[j].qual); free(w2[j].info);
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

/***********************************
 * High-level error correction API *
 ***********************************/

#define DEFAULT_QUAL 20

int fm6_api_correct(int kmer, int64_t l, char *_seq, char *_qual)
{
	char *qual;
	int64_t i, cnt[2];
	int j, *info;
	rld_t *e;
	fmecopt_t opt;
	shash_t **solid;
	char **seq2, **qual2;
	fmintv_t *top;

	// set correction parameters
	opt.w = kmer > 0? kmer : 19;
	opt.min_occ = 3;
	opt.keep_bad = 1; opt.is_paired = 0;
	opt.max_corr = 0.3;
	compute_SUF(opt.w > 15? opt.w - 15 : 1);
	// build FM-index; initialize the k-mer hash table
	assert(_seq[l-1] == 0); // must be NULL terminated
	e = fm6_build2(l, _seq);
	qual = _qual? _qual : malloc(l);
	if (_qual == 0)
		for (i = 0; i < l; ++i)
			qual[i] = DEFAULT_QUAL + 33;
	solid = malloc(SUF_NUM * sizeof(void*));
	for (i = 0; i < SUF_NUM; ++i) solid[i] = kh_init(solid);
	// collect solid k-mers
	top = fm6_traverse(e, SUF_LEN);
	for (i = 0; i < SUF_NUM; ++i)
		ec_collect(e, &opt, SUF_LEN, &top[i], solid[i], cnt);
	free(top);
	// correct errors
	seq2  = malloc(sizeof(void*) * e->mcnt[1] / 2); // NB: e->mcnt[1] equals twice of the number of sequences in _seq
	qual2 = malloc(sizeof(void*) * e->mcnt[1] / 2);
	info  = calloc(e->mcnt[1] / 2, sizeof(int));
	seq2[0] = _seq; qual2[0] = qual;
	for (i = j = 1; i < l - 1; ++i)
		if (_seq[i] == 0)
			seq2[j] = &_seq[i + 1], qual2[j++] = &qual[i + 1];
	ec_fix(e, &opt, solid, e->mcnt[1]/2, seq2, qual2, info);
	free(seq2); free(qual2); free(info);
	// free
	for (i = 0; i < SUF_NUM; ++i) kh_destroy(solid, solid[i]);
	free(solid);
	rld_destroy(e);
	if (_qual == 0) free(qual);
	return 0;
}
