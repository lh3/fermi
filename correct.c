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

typedef struct {
	uint32_t key;
	uint16_t b1:2, b2:2, ratio:6, d2:3, d3:3;
} __attribute__ ((__packed__)) solid1_t;

#define solid_hash(a) (__ac_Wang_hash((a).key))
#define solid_eq(a, b) ((a).key == (b).key)
#define solid_set_key(p, x) ((p)->key = (x))
#include "khash.h"
KHASH_INIT(solid, solid1_t, char, 0, solid_hash, solid_eq)
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
	int i, shift = (opt->w - len - 1) * 2;
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
			solid1_t zz;
			int max_c, max_c2, absent;
			uint64_t max, max2, rest, key, sum;
			double r;
			for (c = 1, max = max2 = 0, max_c = max_c2 = 6, sum = 0; c <= 4; ++c) {
				if (ok[c].x[2] > max) max2 = max, max_c2 = max_c, max = ok[c].x[2], max_c = c;
				else if (ok[c].x[2] > max2) max2 = ok[c].x[2], max_c2 = c;
				sum += ok[c].x[2];
			}
			if (max < opt->min_occ) continue; // then in the following max_c<6
			++cnt[0];
			rest = sum - max;
			r = rest == 0? max : (double)max / rest;
			if (r > 63.) r = 63.; // we have maximally 6 bits of information (i.e. [0,63])
			if (rest <= 7 && r >= opt->min_occ) ++cnt[1];
			for (i = 0, key = 0; i < str.l; ++i)
				key = (uint64_t)str.s[i]<<shift | key>>2;
			solid_set_key(&zz, key);
			zz.b1 = max_c - 1; zz.b2 = max_c2 - 1;
			zz.d2 = rest < 7? rest : 7;
			zz.d3 = sum - max - max2 < 7? sum - max - max2 : 7;
			zz.ratio = (int)(r + .499);
			kh_put(solid, solid, zz, &absent);
			assert(absent);
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
	ku64_v stack;
} fixaux_t;

#define F_NOHIT     0x1
#define F_BEST      0x2
#define F_CHECKED   0x4
#define F_CORRECTED 0x8

#define MAX_QUAL    41

static inline void save_state(fixaux_t *fa, const ku128_t *p, int c, int score, int shift, int flag, int qual)
{
	ku128_t w;
	uint64_t *q;
	int pos = (p->y&0xfffff) - 1;
	assert(pos < 0x100000);
	assert(c >= 0 && c < 4);
	qual = qual < MAX_QUAL? qual : MAX_QUAL;
	score = score > 0? score : 0;
	// update heap; the structure of w.y -- score:16, pos_in_stack:28, seq_pos:20
	w.x = (uint64_t)c<<shift | p->x>>2;
	w.y = (uint64_t)((p->y>>48) + score)<<48 | fa->stack.n<<20 | pos;
	kv_push(ku128_t, fa->heap, w);
	ks_heapup_128y(fa->heap.n, fa->heap.a);
	// update stack; seq_pos:20, qual:8, flag:4, base:4, parent_pos_in_stack:28
	kv_pushp(uint64_t, fa->stack, &q);
	*q = (uint64_t)pos<<44 | (uint64_t)qual<<36 | (uint64_t)flag<<32 | (uint32_t)c<<28 | (p->y>>20&0xfffffff);
}

#define RATIO_FACTOR  10
#define DIFF_FACTOR   13
#define MAX_HEAP      256
#define MAX_SC_DIFF   60
#define MISS_PENALTY  10
#define MIN_OCC       5
#define MIN_OCC_RATIO 0.8

typedef struct {
	int n_hash_hit;
	int n_hash_qry;
	int l_hash_cov;
	int min_penalty;
	int pen_diff;
	int l_cov;
} ecstat1_t;

static inline void ec_cal_penalty(const solid1_t *p, int penalty[2], int base)
{
	int tmp, d1 = p->d2? p->d2 * p->ratio : p->ratio;
	penalty[0] = ((d1 < 7? d1 : 7) - p->d2) * DIFF_FACTOR;
	tmp = p->d2? p->ratio * RATIO_FACTOR : 10000;
	penalty[0] = tmp < penalty[0]? tmp : penalty[0];
	penalty[0] = penalty[0] > 1? penalty[0] : 1;
	if (p->b1 != base && p->b2 != base) {
		d1 += p->d2 - p->d3;
		penalty[1] = ((d1 < 7? d1 : 7) - p->d3) * DIFF_FACTOR;
		tmp = p->d3? (int)((double)d1 / p->d3 * RATIO_FACTOR + .499) : 10000;
		penalty[1] = tmp < penalty[1]? tmp : penalty[1];
		penalty[1] = penalty[1] > 1? penalty[1] : 1;
		//penalty[0] += penalty[1];
	} else penalty[1] = 0;
	penalty[0] = penalty[0] < MAX_QUAL? penalty[0] : MAX_QUAL;
	penalty[1] = penalty[1] < MAX_QUAL? penalty[1] : MAX_QUAL;
}

typedef struct {
	uint8_t ob, oq, cb, cq;
} ecbase_t;

typedef kvec_t(ecbase_t) ecseq_t;

static void ec_seq_gen(ecseq_t *s, int l_seq, char *seq, char *qual)
{
	int i;
	kv_resize(ecbase_t, *s, l_seq);
	s->n = l_seq;
	for (i = 0; i < l_seq; ++i) {
		ecbase_t *b = &s->a[i];
		b->ob = b->cb = seq_nt6_table[(int)seq[i]];
		b->oq = qual[i] - 33;
		b->cq = 0;
	}
}

static void ec_seq_rev(ecseq_t *s)
{
	int i;
	for (i = 0; i < s->n>>1; ++i) {
		ecbase_t tmp;
		tmp = s->a[i]; s->a[i] = s->a[s->n - 1 - i]; s->a[s->n - 1 - i] = tmp;
	}
	for (i = 0; i < s->n; ++i) {
		s->a[i].ob = fm6_comp(s->a[i].ob);
		s->a[i].cb = fm6_comp(s->a[i].cb);
	}
}

static void ec_fix1(const fmecopt_t *opt, shash_t *const* solid, ecseq_t *s, fixaux_t *fa, ecstat1_t *es)
{
	int i, l, shift = (opt->w - 1) << 1, n_rst = 0;
	ku128_t z, rst[2];

	memset(es, 0, sizeof(ecstat1_t));
	if (s->n <= opt->w) return;
	fa->heap.n = fa->stack.n = 0;
	kv_push(uint64_t, fa->stack, 0);
	// get the initial k-mer
	for (i = s->n - 1, l = 0, z.x = 0; i > 0 && l < opt->w; --i)
		if (s->a[i].cb == 5) z.x = 0, l = 0;
		else z.x = (uint64_t)(s->a[i].cb-1)<<shift | z.x>>2, ++l;
	if (i == 0) return; // no good k-mer
	// the first element in the heap
	z.y = i + 1;
	kv_push(ku128_t, fa->heap, z);
	// traverse
	while (fa->heap.n) {
		const shash_t *h;
		khint_t k;
		solid1_t zz;
		int cq, qual, penalty[2];
		// get the best so far
		z = fa->heap.a[0]; // $z is the best
		fa->heap.a[0] = kv_pop(fa->heap);
		ks_heapdown_128y(0, fa->heap.n, fa->heap.a);
		if ((z.y&0xfffff) == 0) {
			rst[n_rst++] = z;
			if (n_rst == 2) break;
			continue;
		}
		if (n_rst && (int)(z.y>>48) > (int)(rst[0].y>>48) + MAX_SC_DIFF) break;
		i = (z.y&0xfffff) - 1; // $i points to the base to be checked
		qual = s->a[i].oq < MAX_QUAL? s->a[i].oq : MAX_QUAL;
		if (qual < 3) qual = 3;
		if (fm_verbose >= 5) fprintf(stderr, "pop\tsc=%d\ti=%d\t%c%d\n", (int)(z.y>>48), i, "$ACGTN"[(int)s->a[i].cb], qual);
		// check the hash table
		h = solid[z.x & (SUF_NUM - 1)];
		solid_set_key(&zz, z.x >> (SUF_LEN<<1));
		k = kh_get(solid, h, zz);
		++es->n_hash_qry;
		if (k != kh_end(h)) { // this (k+1)-mer has more than opt->min_occ occurrences
			++es->n_hash_hit;
			if (s->a[i].cb != kh_key(h, k).b1 + 1) { // the read base is different from the best base
				solid1_t *p = &kh_key(h, k);
				ec_cal_penalty(p, penalty, s->a[i].cb - 1);
				// if we have too many possibilities, keep the better path among the two
				if (s->a[i].cb != 5 && (fa->heap.n + 2 <= MAX_HEAP || penalty[0] < qual)) { // the read path
					cq = qual > penalty[0]? qual - penalty[0] : 2;
					save_state(fa, &z, s->a[i].cb - 1, penalty[0], shift, F_CHECKED, cq);
				}
				if (s->a[i].cb == 5 || fa->heap.n + 2 <= MAX_HEAP || penalty[0] > qual) { // the stack path
					cq = penalty[0] > qual? penalty[0] - qual : 2;
					save_state(fa, &z, p->b1, qual, shift, F_CHECKED|F_CORRECTED, cq);
				}
				if (fm_verbose >= 5) fprintf(stderr, "cmp\ti=%d\t%c%d => %c%d?\n", i, "$ACGTN"[(int)s->a[i].cb], qual, "ACGT"[p->b1], penalty[0]);
			} else { // the read base is the same as the best base
				ku128_t z0 = z;
				solid1_t *p = &kh_key(h, k);
				int i0 = i;
				if (p->d2 <= 0 && opt->step > 1) {
					/* The following is a complex block. In principle, we can
					 * skip it altogether, which will yield a little better
					 * accuracy. However, without this block, we will need
					 * much more hash table lookups, which is the bottleneck.
					 */
					int last_penalty, occ_last = p->d2? p->d2 * (p->ratio+1) : p->ratio;
					ec_cal_penalty(p, penalty, s->a[i0].cb - 1);
					last_penalty = penalty[0];
					while (i0 > 0) {
						for (i = (z.y&0xfffff) - 1, l = 0; i >= 1 && l < opt->step && s->a[i-1].cb < 5 && s->a[i-1].cb == s->a[i-1].ob; --i, ++l)
							z.x = (uint64_t)(s->a[i].cb-1)<<shift | z.x>>2; // look opt->step mer ahead
						if (l <= 1) break; // making no progress or only one step further; stop
						h = solid[z.x & (SUF_NUM - 1)];
						solid_set_key(&zz, z.x >> (SUF_LEN<<1));
						k = kh_get(solid, h, zz);
						++es->n_hash_qry;
						es->n_hash_hit += (k != kh_end(h));
						if (k != kh_end(h) && s->a[i].cb == kh_key(h, k).b1 + 1) { // in the hash table and the read base is the best
							solid1_t *p = &kh_key(h, k);
							int pen_save, parent, occ = p->d2? p->d2 * (p->ratio+1) : p->ratio;
							if (fm_verbose >= 5) fprintf(stderr, "jump\ti=%d\t%c%d\t%d\n", i, "$ACGTN"[(int)s->a[i].cb], s->a[i].oq, occ);
							if (p->d2 > 1 || occ < MIN_OCC || occ < occ_last * MIN_OCC_RATIO) break; // if occ is not good enough; stop
							// update stack
							ec_cal_penalty(p, penalty, s->a[i].cb - 1);
							pen_save = last_penalty < penalty[0]? last_penalty : penalty[0];
							last_penalty = penalty[0];
							parent = z.y>>20 & 0xfffffff;
							while (--l >= 0) {
								uint64_t *q;
								int pos = (z.y&0xfffff) - 1 - l;
								int qual = pen_save > s->a[pos].oq? pen_save : s->a[pos].oq;
								kv_pushp(uint64_t, fa->stack, &q);
								*q = (uint64_t)pos<<44 | (uint64_t)qual<<36 | (uint64_t)F_BEST<<32 | (uint32_t)(s->a[pos].cb-1)<<28 | parent;
								parent = fa->stack.n - 1;
							}
							// prepare for the next round of iteration
							z.y = z.y>>48<<48 | (uint64_t)parent<<20 | (i + 1);
							z0 = z; i0 = i;
							occ_last = occ;
						} else break;
					}
				}
				ec_cal_penalty(p, penalty, s->a[i0].cb - 1);
				save_state(fa, &z0, s->a[i0].cb - 1, 0, shift, F_BEST, penalty[0]);
			}
		} else save_state(fa, &z, s->a[i].cb - 1, MISS_PENALTY + (MAX_QUAL - qual), shift, F_NOHIT, 0);
	}
	assert(n_rst == 1 || n_rst == 2);
	es->min_penalty = rst[0].y >> 48;
	es->pen_diff = n_rst == 1? MAX_SC_DIFF : (int)(rst[1].y>>48) - (int)(rst[0].y>>48);
	assert(es->pen_diff >= 0);
	es->pen_diff = es->pen_diff < MAX_SC_DIFF? es->pen_diff : MAX_SC_DIFF;
	// backtrack
	l = rst[0].y>>20&0xfffffff;
	while (l) {
		i = fa->stack.a[l]>>44;
		s->a[i].cb = ((uint32_t)fa->stack.a[l]>>28) + 1;
		s->a[i].cq = fa->stack.a[l]>>36;
		l = fa->stack.a[l] & 0xfffffff; // take the lowest 28 bits, the parent position in the stack
	}
}

typedef struct {
	uint32_t q_corr:8, pen_diff:8, p_cov:7, no_hit:1, p_corr:7, dummy:1;
} ecinfo1_t;

static uint64_t ec_fix(const rld_t *e, const fmecopt_t *opt, shash_t *const* solid, int n_seqs, char **seq, char **qual, ecinfo1_t *info)
{
	int i;
	uint64_t n_query = 0;
	ecseq_t s;
	fixaux_t fa;
	fmintv_v prev, curr;
	fmintv6_v tmp;
	fmec2opt_t opt2;

	kv_init(s);
	kv_init(prev); kv_init(curr); kv_init(tmp);
	fmc_opt_init(&opt2);
	memset(&fa, 0, sizeof(fixaux_t));
	for (i = 0; i < n_seqs; ++i) {

		if (1) {
			fm6_ec2_core(&opt2, e, strlen(seq[i]), seq[i], qual[i], &prev, &curr, &tmp);
			continue;
		}

		char *si = seq[i], *qi = qual[i];
		ecinfo1_t *ii = &info[i];
		ecstat1_t es[2];

		memset(es, 0, sizeof(ecstat1_t) * 2);
		ec_seq_gen(&s, strlen(si), si, qi);
		ec_seq_rev(&s);
		if (fm_verbose >= 5) fprintf(stderr, "=== index %d, forward ===\n", i);
		ec_fix1(opt, solid, &s, &fa, &es[0]); // 0x7fff0000 if no correction; 0xffff if too short
		n_query += es[0].n_hash_qry;
		ii->pen_diff = es[0].pen_diff;
		ec_seq_rev(&s);
		if (es[0].n_hash_qry) { // then we need to correct in the reverse direction
			if (fm_verbose >= 5) fprintf(stderr, "=== index %d, reverse ===\n", i);
			ec_fix1(opt, solid, &s, &fa, &es[1]);
			n_query += es[1].n_hash_qry;
			ii->pen_diff = ii->pen_diff < es[1].pen_diff? ii->pen_diff : es[1].pen_diff;
		}
		if (es[0].n_hash_hit || es[1].n_hash_hit) {
			int j, n_corr = 0, q_corr = 0, l_cov = 0;
			for (j = 0; j < s.n; ++j) {
				ecbase_t *p = &s.a[j];
				int is_diff = (p->ob != p->cb);
				si[j] = is_diff? "$acgtn"[(int)p->cb] : "$ACGTN"[(int)p->cb];
				if (p->cq) ++l_cov, qi[j] = (p->cq>>1<<1|1) + 33;
				else qi[j] = (p->oq>>1<<1) + 33;
				n_corr += is_diff;
				q_corr += is_diff? p->oq : 0;
			}
			ii->no_hit = 0;
			ii->q_corr = q_corr < 255? q_corr : 255;
			ii->p_cov = (int)(100. * l_cov / s.n + .499);
			ii->p_corr = (int)(100. * n_corr / s.n + .499);
		} else ii->no_hit = 1;
	}
	free(s.a);
	free(fa.heap.a); free(fa.stack.a);
	free(prev.a); free(curr.a); free(tmp.a);
	return n_query;
}

/***************************
 * Dump/restore K-mer hash *
 ***************************/

static void dump_hash(const char *fn, const fmecopt_t *opt, shash_t **solid)
{
	gzFile fp;
	int i;
	fp = strcmp(fn, "-")? gzopen(fn, "w1") : gzdopen(fileno(stdout), "w1");
	gzwrite(fp, opt, sizeof(fmecopt_t));
	gzwrite(fp, &SUF_LEN, sizeof(int));
	for (i = 0; i < SUF_NUM; ++i) {
		gzwrite(fp, solid[i], sizeof(shash_t)); // we don't actually need to dump pointers, but it does not matter too much
		gzwrite(fp, solid[i]->flags, __ac_fsize(solid[i]->n_buckets) * 4);
		gzwrite(fp, solid[i]->keys, solid[i]->n_buckets * sizeof(solid1_t));
	}
	gzclose(fp);
}

static shash_t **restore_hash(const char *fn, fmecopt_t *opt)
{
	gzFile fp;
	int i;
	shash_t **solid;
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	gzread(fp, opt, sizeof(fmecopt_t));
	gzread(fp, &SUF_LEN, sizeof(int));
	compute_SUF(SUF_LEN);
	solid = calloc(SUF_NUM, sizeof(void*));
	for (i = 0; i < SUF_NUM; ++i) {
		solid[i] = malloc(sizeof(shash_t));
		gzread(fp, solid[i], sizeof(shash_t));
		solid[i]->flags = malloc(__ac_fsize(solid[i]->n_buckets) * 4);
		gzread(fp, solid[i]->flags, __ac_fsize(solid[i]->n_buckets) * 4);
		solid[i]->keys = malloc(solid[i]->n_buckets * sizeof(solid1_t));
		gzread(fp, solid[i]->keys, solid[i]->n_buckets * sizeof(solid1_t));
	}
	gzclose(fp);
	return solid;
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
	int n_seqs;
	ecinfo1_t *info;
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

int fm6_ec_correct(const rld_t *e, fmecopt_t *opt, const char *fn, int _n_threads, const char *fn_hash)
{
	int j, n_threads;
	int64_t i, cnt[2];
	shash_t **solid;
	pthread_t *tid;

	if (!fn_hash) {
		if (opt->w < 0) { // determine k-mer
			opt->w = (int)(log(e->mcnt[0]) / log(4) + 8.499);
			if (opt->w >= MAX_KMER) opt->w = MAX_KMER;
			if (fm_verbose >= 3)
				fprintf(stderr, "[M::%s] set k-mer length to %d\n", __func__, opt->w);
		}
		compute_SUF(opt->w > 16? opt->w - 16 : 1);
		solid = calloc(SUF_NUM, sizeof(void*));
		for (j = 0; j < SUF_NUM; ++j) solid[j] = kh_init(solid);
	} else {
		fmecopt_t old = *opt;
		solid = restore_hash(fn_hash, opt);
		opt->step = old.step; opt->trim_l = old.trim_l;
	}

	assert(_n_threads <= SUF_NUM);
	tid = (pthread_t*)calloc(_n_threads, sizeof(pthread_t));
	cnt[0] = cnt[1] = 0;

	if (!fn_hash) { // initialize and launch worker1
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
		for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], 0, worker1, w1 + j);
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

	if (!fn) { // no sequence file is given; dump the hash table only without correction
		dump_hash("-", opt, solid);
	} else { // initialize and launch worker2
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
			w2[j].info = calloc(max_seqs, sizeof(ecinfo1_t));
		}
		for (;;) {
			int ret;
			worker2_t *w;
			ret = kseq_read(seq);
			if (ret < 0 || (id && id%BATCH_SIZE == 0)) {
				uint64_t n_query = 0, n_seq = 0;
				for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], 0, worker2, w2 + j);
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
					int tmp = 0;
					ecinfo1_t *info;
					w = &w2[k%n_threads];
					info = &w->info[w->n_seqs];
					out.l = 0;
					kputc('@', &out); kputl(k, &out);
					kputc('_', &out); kputw(info->no_hit, &out);
					kputc('_', &out); kputw(info->pen_diff, &out);
					kputc('_', &out); kputw(info->q_corr, &out);
					kputc('_', &out); kputw(info->p_corr, &out);
					kputc('_', &out); kputw(info->p_cov, &out);
					kputc('\n', &out);
					tmp = strlen(w->seq[w->n_seqs]);
					if (opt->trim_l && opt->trim_l < tmp) tmp = opt->trim_l;
					kputsn(w->seq[w->n_seqs], tmp, &out);
					kputsn("\n+\n", 3, &out); kputsn(w->qual[w->n_seqs], tmp, &out);
					puts(out.s);
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

/* Failing examples:

   AATGTAGGGATTGTTTTAATTCCCGCTCCCTTATGGGAGAAGGTTAACGCTCGGGTAACCCTTGCCGAATGTAGGCCGGATAAGGCGTTTACGCtGCATCC (NC_000913)
   IHIIIHHIIIIIFIIIIIIIIIIIIIHIIIIIIIIIIGBIIGGHIHHIIHGIIFHEIIIIHIIHIIGBGIIIC<@BEEGB3EC<DED@CEDDDB43C?AAC

   The lowercase 't' is an error but not corrected as the 22mer preceding it is a repeat.
 */

/***********************************
 * High-level error correction API *
 ***********************************/

#define DEFAULT_QUAL 20

int fm6_api_correct(int kmer, int64_t l, char *_seq, char *_qual)
{
	char *qual;
	int64_t i, cnt[2];
	int j;
	ecinfo1_t *info;
	rld_t *e;
	fmecopt_t opt;
	shash_t **solid;
	char **seq2, **qual2;
	fmintv_t *top;

	// set correction parameters
	opt.w = kmer > 0? kmer : 19;
	opt.min_occ = 3;
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
	info  = calloc(e->mcnt[1] / 2, sizeof(ecinfo1_t));
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
