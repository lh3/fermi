#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <stdio.h>
#include <zlib.h>
#include "khash.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#ifndef _BFC_EXT_TABLE
static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};
#else
extern unsigned char seq_nt6_table[128];
#endif

static inline uint32_t hash_32(uint32_t key)
{
	key += ~(key << 15);
	key ^=  (key >> 10);
	key +=  (key << 3);
	key ^=  (key >> 6);
	key += ~(key << 11);
	key ^=  (key >> 16);
	return key;
}

typedef struct {
	uint32_t x;
	uint8_t y;
} __attribute__ ((__packed__)) kmer1_t;

#define km_hash1(a) (hash_32((a).x))
#define km_eq1(a, b) ((a).x == (b).x)
KHASH_INIT(km1, kmer1_t, char, 0, km_hash1, km_eq1)

#define km_hash2(a) (hash_32((a).x>>2))
#define km_eq2(a, b) ((a).x>>2 == (b).x>>2)
KHASH_INIT(km2, kmer1_t, char, 0, km_hash2, km_eq2)
typedef khash_t(km2) shash_t;

static inline uint64_t bfc_hash(uint64_t key, uint64_t mask)
{
	key = (key + ~(key << 32)) & mask;
	key = (key ^ (key >> 22)) & mask;
	key = (key + (key << 13)) & mask;
	key = (key ^ (key >> 24)) & mask;
	key = (key + ~(key << 27)) & mask;
	key = (key ^ (key >> 31)) & mask;
	return key;
}

static inline uint64_t bfc_hash_inv(uint64_t key, uint64_t mask)
{ // for inversion, see also: http://naml.us/blog/tag/invertible
	uint64_t tmp;

	tmp = (key ^ key >> 31) & mask;
	key = (key ^ tmp >> 31) & mask;

	tmp = (key - (~key << 27)) & mask;
	key = (key - (~tmp << 27) + 1) & mask;

	tmp = (key ^ key >> 24) & mask;
	key = (key ^ tmp >> 24) & mask;

	tmp = (key - (key << 13)) & mask;
	tmp = (key - (tmp << 13)) & mask;
	tmp = (key - (tmp << 13)) & mask;
	key = (key - (tmp << 13)) & mask;

	tmp = (key ^ key >> 22) & mask;
	tmp = (key ^ tmp >> 22) & mask;
	key = (key ^ tmp >> 22) & mask;

	tmp = (key - (~key << 32)) & mask;
	key = (key - (~tmp << 32) + 1) & mask;

	return key;
}

static inline uint64_t bfc_hash2(uint64_t key)
{
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

/***********************
 *** Large bit array ***
 ***********************/

#define BF_SHIFT 20
#define BF_MASK  ((1<<BF_SHIFT) - 1)

typedef struct {
	int n, n_bits;
	uint64_t **a;
} barray_t;

static inline barray_t *ba_init(int n_bits)
{
	int i;
	barray_t *ba;
	ba = calloc(1, sizeof(barray_t));
	ba->n_bits = n_bits;
	ba->n = ((1<<n_bits) + (1<<BF_SHIFT) - 1) >> BF_SHIFT;
	ba->a = calloc(ba->n, sizeof(void*));
	for (i = 0; i < ba->n; ++i)
		ba->a[i] = calloc(1<<BF_SHIFT>>6, 8);
	return ba;
}

static inline void ba_destroy(barray_t *ba)
{
	int i;
	for (i = 0; i < ba->n; ++i) free(ba->a[i]);
	free(ba->a); free(ba);
}

static inline int ba_set(barray_t *ba, uint64_t x)
{
	uint64_t *p = &ba->a[x>>BF_SHIFT][(x&BF_MASK)>>6];
	uint64_t r, z = 1ULL<<(x&63);
	r = __sync_fetch_and_or(p, z);
	return r>>(x&63)&1;
}

static inline int ba_get(const barray_t *ba, uint64_t x)
{
	return ba->a[x>>BF_SHIFT][(x&BF_MASK)>>6]>>(x&63)&1;
}

/************************
 *** K-mer hash table ***
 ************************/

typedef struct {
	int k;
	barray_t *a;
	void **h;
} bfc_hash_t;

bfc_hash_t *bfc_hash_init(int k, int ba_bits)
{
	bfc_hash_t *b;
	int nh, i;
	b = calloc(1, sizeof(bfc_hash_t));
	b->k = k;
	b->a = ba_init(ba_bits);
	nh = k > 32? 1<<(k-32) : 1;
	b->h = calloc(nh, sizeof(void*));
	for (i = 0; i < nh; ++i)
		b->h[i] = kh_init(km1);
	return b;
}

void bfc_hash_destroy(bfc_hash_t *b)
{
	int i, nh;
	nh = b->k > 32? 1<<(b->k-32) : 1;
	for (i = 0; i < nh; ++i)
		kh_destroy(km1, (khash_t(km1)*)b->h[i]);
	free(b->h);
	ba_destroy(b->a);
	free(b);
}

uint64_t bfc_hash_count1(const bfc_hash_t *b)
{
	int i, nh = b->k > 32? 1<<(b->k-32) : 1;
	uint64_t size = 0;
	for (i = 0; i < nh; ++i)
		size += kh_size((khash_t(km1)*)b->h[i]);
	return size;
}

void bfc_hash_put(bfc_hash_t *b, uint64_t x)
{
	uint64_t bmask = (1ULL << b->a->n_bits) - 1;
	khash_t(km1) *h;
	int r = 0;
	kmer1_t tmp;
	khint_t k;

	r += ba_set(b->a, bfc_hash(x, (1ULL<<b->k) - 1) & bmask);
	r += ba_set(b->a, bfc_hash2(x) & bmask);
	r += ba_set(b->a, bfc_hash_inv(x, (1ULL<<b->k) - 1) & bmask);
	if (r < 3) return;
	tmp.x = (uint32_t)x; tmp.y = 0;
	h = (khash_t(km1)*)b->h[x>>32];
	while (__sync_lock_test_and_set(&h->lock, 1)) while (h->lock);
	k = kh_put(km1, h, tmp, &r);
	if (kh_key(h, k).y < 255) ++kh_key(h, k).y;
	__sync_lock_release(&h->lock);
}

/**********************
 *** Collect K-mers ***
 **********************/

typedef struct {
	int min_q, k, ba_bits;
	int n_threads;
	int chunk_size;
	int def_qual;
} bfc_opt_t;

typedef struct {
	int l;
	uint8_t *s, *q;
} seq1_t;

typedef struct {
	int n, start, step;
	seq1_t *seq;
	bfc_hash_t *h;
	const bfc_opt_t *opt;
} worker1_t;

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

static void add_seq1(const bfc_opt_t *opt, bfc_hash_t *h, int l, const uint8_t *s, const uint8_t *q)
{
	int k, i, shift = opt->k << 1;
	uint64_t f, r, mask = (1ULL << ((opt->k+1)<<1)) - 1;
	for (f = r = 0, i = k = 0; i < l; ++i) {
		if (s[i] < 4) {
			f = (f<<2 | s[i]) & mask; // forward
			r = r>>2 | (uint64_t)(3-s[i])<<shift; // reverse
			if (++k >= opt->k+1 && q[i] >= opt->min_q) {
				bfc_hash_put(h, f);
				bfc_hash_put(h, r);
			}
		} else k = 0;
	}
}

static void *worker1(void *data)
{
	worker1_t *w = (worker1_t*)data;
	int i, j;
	for (i = w->start; i < w->n; i += w->step) {
		seq1_t *s = &w->seq[i];
		if (s->q == 0) {
			s->q = malloc(s->l);
			memset(s->q, 33 + w->opt->def_qual, s->l);
		}
		for (j = 0; j < s->l; ++j) {
			s->s[j] = s->s[j] > 127? 4 : seq_nt6_table[s->s[j]] - 1;
			s->q[j] -= 33;
		}
		add_seq1(w->opt, w->h, s->l, s->s, s->q);
		free(s->s); free(s->q);
	}
	return 0;
}

static void collect_core(const bfc_opt_t *opt, bfc_hash_t *h, int n, seq1_t *seq)
{
	worker1_t *w;
	pthread_t *tid;
	int j;
	tid = malloc(opt->n_threads * sizeof(pthread_t));
	w = calloc(opt->n_threads, sizeof(worker1_t));
	for (j = 0; j < opt->n_threads; ++j) {
		w[j].n = n; w[j].start = j; w[j].step = opt->n_threads;
		w[j].seq = seq; w[j].h = h; w[j].opt = opt;
	}
	for (j = 0; j < opt->n_threads; ++j) pthread_create(&tid[j], 0, worker1, &w[j]);
	for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
	free(w); free(tid);
}

void bfc_opt_init(bfc_opt_t *opt)
{
	opt->min_q = 17;
	opt->chunk_size = 1<<26;
	opt->n_threads = 1;
	opt->k = 23;
	opt->ba_bits = 23;
	opt->def_qual = 20;
}

bfc_hash_t *bfc_collect(const bfc_opt_t *opt, const char *fn)
{
	void *ko;
	int fd, n, m, tot_len;
	gzFile fp;
	kseq_t *ks;
	seq1_t *seq;
	bfc_hash_t *h = 0;

	ko = kopen(fn, &fd);
	fp = gzdopen(fd, "r");
	ks = kseq_init(fp);
	h = bfc_hash_init((opt->k + 1) * 2, opt->ba_bits);
	n = m = tot_len = 0;
	seq = 0;
	while (kseq_read(ks) >= 0) {
		if (tot_len >= opt->chunk_size) {
			collect_core(opt, h, n, seq);
			tot_len = n = 0;
		}
		if (n == m) {
			m = m? m<<1 : 1024;
			seq = realloc(seq, m * sizeof(seq1_t));
		}
		seq[n].l = ks->seq.l;
		seq[n].s = malloc(ks->seq.l);
		memcpy(seq[n].s, ks->seq.s, ks->seq.l);
		if (ks->qual.l) {
			seq[n].q = malloc(ks->seq.l);
			memcpy(seq[n].q, ks->qual.s, ks->seq.l);
		} else seq[n].q = 0;
		++n; tot_len += ks->seq.l;
	}
	collect_core(opt, h, n, seq);
	free(seq);
	kseq_destroy(ks);
	gzclose(fp);
	kclose(ko);
	fprintf(stderr, "[M::%s] collected approximately %ld non-unique k-mers\n", __func__, (long)bfc_hash_count1(h));
	return h;
}

int main(int argc, char *argv[])
{
	int c;
	bfc_opt_t opt;
	bfc_hash_t *h;

	bfc_opt_init(&opt);
	while ((c = getopt(argc, argv, "k:t:b:")) >= 0) {
		if (c == 'k') opt.k = atoi(optarg);
		else if (c == 't') opt.n_threads = atoi(optarg);
		else if (c == 'b') opt.ba_bits = atoi(optarg);
	}

	if (optind + 1 > argc) {
		return 1;
	}

	h = bfc_collect(&opt, argv[optind]);
	bfc_hash_destroy(h);
	return 0;
}
