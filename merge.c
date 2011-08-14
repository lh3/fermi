#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "fermi.h"
#include "rld.h"

double cputime();
double rssmem();

typedef struct {
	rlditr_t itr;
	int c;
	int64_t l;
} rlditr2_t;

// This is a clever version of rld_enc(): adjacent runs are guaranteed to be different.
inline void rld_enc2(rld_t *e, rlditr2_t *itr2, int64_t l, int c)
{
	if (l == 0) return;
	if (itr2->c != c) {
		if (itr2->l) rld_enc(e, &itr2->itr, itr2->l, itr2->c);
		itr2->l = l; itr2->c = c;
	} else itr2->l += l;
}
// take k symbols from e0 and write it to e; l0 number of c0 symbols are pending before writing
static inline void dec_enc(rld_t *e, rlditr2_t *itr, const rld_t *e0, rlditr_t *itr0, int64_t *l0, int *c0, int64_t k)
{
	if (*l0 >= k) { // there are more pending symbols
		rld_enc2(e, itr, k, *c0);
		*l0 -= k; // l0-k symbols remains
	} else { // use up all pending symbols
		int c = -1; // to please gcc
		int64_t l;
		rld_enc2(e, itr, *l0, *c0); // write all pending symbols
		k -= *l0;
		for (; k > 0; k -= l) { // we always go into this loop because l0<k
			l = rld_dec(e0, itr0, &c, 1);
			assert(l); // the e0 stream should not be finished
			rld_enc2(e, itr, k < l? k : l, c);
		}
		*l0 = -k; *c0 = c;
	}
}

#include "khash.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint64_t)
#define h64_eq(a, b) ((a)>>32 == (b)>>32)
#define h64_hash(a) ((a)>>32)
KHASH_INIT(h64, uint64_t, char, 0, h64_hash, h64_eq)

#define BLOCK_BITS 16
#define BLOCK_MASK ((1u<<BLOCK_BITS) - 1)
#define BLOCK_SHIFT (64 - BLOCK_BITS)
#define BLOCK_CMASK ((1ll<<BLOCK_SHIFT) - 1)

#define MSG_SIZE 10000000

typedef struct {
	int n;
	khash_t(h64) **h;
} gaphash_t;

static inline gaphash_t *init_gaphash(uint64_t n)
{
	int i;
	gaphash_t *h;
	h = calloc(1, sizeof(gaphash_t));
	h->n = (n + BLOCK_MASK) >> BLOCK_BITS;
	h->h = malloc(h->n * sizeof(void*));
	for (i = 0; i < h->n; ++i)
		h->h[i] = kh_init(h64);
	return h;
}

static inline void insert_to_hash(gaphash_t *h, uint64_t j)
{
	khint_t k;
	int ret;
	khash_t(h64) *g = h->h[j>>BLOCK_BITS];
	k = kh_put(h64, g, (j&BLOCK_MASK)<<BLOCK_SHIFT|1, &ret);
	if (ret == 0) ++kh_key(g, k); // when the key is present, the key in the hash table will not be overwritten by put()
}

static void *compute_gap(const rld_t *e0, const rld_t *e1, int use_hash)
{
	uint64_t k, *ok, i, x, *bits = 0, n_processed = 1;
	int c = 0;
	double t = cputime();
	gaphash_t *h = 0;

	ok = alloca(8 * e0->asize);
	x = e1->mcnt[1];
	k = --x; // get the last sentinel of e1
	i = e0->mcnt[1] - 1; // to modify gap[j]
	if (use_hash) {
		h = init_gaphash(e0->mcnt[0]);
		insert_to_hash(h, i);
	} else {
		uint64_t n_bits = e0->mcnt[0] + e1->mcnt[0];
		bits = calloc((n_bits + 63) / 64, 8); // we could allocate bits in several blocks to avoid allocating a huge array
		bits[(k+i+1)>>6] |= 1ull<<((k+i+1)&0x3f);
	}
	for (;;) {
		c = rld_rank1a(e1, k, ok);
		if (c == 0) { // sentinel; the order of sentinel has been lost; we have to rely on x to get it back
			k = --x;
			if (x == (uint64_t)-1) break;
			i = e0->mcnt[1] - 1;
		} else {
			k = e1->cnt[c] + ok[c] - 1;
			rld_rank1a(e0, i, ok);
			i = e0->cnt[c] + ok[c] - 1;
		}
		if (use_hash) insert_to_hash(h, i);
		else bits[(k+i+1)>>6] |= 1ull<<((k+i+1)&0x3f);
		if (++n_processed % MSG_SIZE == 0 && fm_verbose >= 3)
			fprintf(stderr, "[M::%s] processed %lld million symbols in %.3f seconds (peak memory: %.3f MB).\n", __func__,
					(long long)n_processed / 1000000, cputime() - t, rssmem());
	}
	return use_hash? (void*)h : (void*)bits;
}

rld_t *fm_merge(rld_t *e0, rld_t *e1, int use_hash)
{
	void *gaparr;
	uint64_t i, n = e0->mcnt[0] + e1->mcnt[0];
	int64_t l0, l1;
	int c0, c1;
	rlditr_t itr0, itr1;
	rlditr2_t itr;
	rld_t *e;

	// compute the gap array
	gaparr = compute_gap(e0, e1, use_hash);
	free(e0->frame); free(e1->frame); // deallocate the rank indexes of e0 and e1; they are not needed any more
	e0->frame = e1->frame = 0;
	// initialize the FM-index to be returned, and all the three iterators
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr.itr, 0);
	itr.l = l0 = l1 = 0; itr.c = c0 = c1 = -1;
	rld_itr_init(e0, &itr0, 0);
	rld_itr_init(e1, &itr1, 0);
	// use the gap array to merge the two BWTs
	if (use_hash) {
		int64_t last = -1;
		int i;
		khint_t k, l;
		gaphash_t *gap = (gaphash_t*)gaparr;
		for (i = 0; i < gap->n; ++i) {
			khash_t(h64) *h = gap->h[i];
			for (l = 0, k = kh_begin(h); k < kh_end(h); ++k)
				if (kh_exist(h, k))
					h->keys[l++] = kh_key(h, k);
			assert(l == kh_size(h));
			free(h->flags);
			h->flags = 0;
			ks_introsort(uint64_t, kh_size(h), h->keys); // actually we do not really need sorting, but this is perhaps faster
			for (k = 0; k < kh_size(h); ++k) {
				uint64_t x = (uint64_t)i<<BLOCK_BITS | h->keys[k]>>BLOCK_SHIFT;
				//printf("gap[%lld]=%lld\n", x, h->keys[k]&BLOCK_CMASK);
				dec_enc(e, &itr, e0, &itr0, &l0, &c0, x - last);
				dec_enc(e, &itr, e1, &itr1, &l1, &c1, h->keys[k]&BLOCK_CMASK);
				last = x;
			}
			kh_destroy(h64, h);
		}
		if (last != e0->mcnt[0] - 1)
			dec_enc(e, &itr, e0, &itr0, &l0, &c0, e0->mcnt[0] - 1 - last);
		free(gap->h);
	} else {
		uint64_t k = 1, *bits = (uint64_t*)gaparr;
		int last = bits[0]&1;
		for (i = 1; i < n; ++i) {
			int c = bits[i>>6]>>(i&0x3f)&1;
			if (c != last) {
				if (last == 0) dec_enc(e, &itr, e0, &itr0, &l0, &c0, k);
				else dec_enc(e, &itr, e1, &itr1, &l1, &c1, k);
				last = c; k = 1;
			} else ++k;
		}
		if (k) {
			if (last == 0) dec_enc(e, &itr, e0, &itr0, &l0, &c0, k);
			else dec_enc(e, &itr, e1, &itr1, &l1, &c1, k);
		}
		free(bits);
	}
	// finalize the merge
	assert(l0 == 0 && l1 == 0); // both e0 and e1 stream should be finished
	rld_enc(e, &itr.itr, itr.l, itr.c); // write the remaining symbols in the iterator
	rld_destroy(e0); rld_destroy(e1);
	rld_enc_finish(e, &itr.itr);
	return e;
}

rld_t *fm_merge_from_SA(rld_t *e0, int len, const uint8_t *T, const int *SA, const int64_t *rank_l)
{
	int64_t l0, last = -1;
	int c0, i;
	rlditr_t itr0;
	rlditr2_t itr;
	rld_t *e;

	free(e0->frame); e0->frame = 0;
	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &itr.itr, 0);
	itr.l = l0 = 0; itr.c = c0 = -1;
	rld_itr_init(e0, &itr0, 0);
	for (i = 0; i < len; ++i) {
		if (rank_l[i] != last) {
			dec_enc(e, &itr, e0, &itr0, &l0, &c0, rank_l[i] - last);
			last = rank_l[i];
		}
		rld_enc2(e, &itr, 1, SA[i]? T[SA[i]-1] : 0);
	}
	if (last != e0->mcnt[0] - 1)
		dec_enc(e, &itr, e0, &itr0, &l0, &c0, e0->mcnt[0] - 1 - last);
	rld_enc(e, &itr.itr, itr.l, itr.c); // write the remaining symbols in the iterator
	rld_destroy(e0);
	rld_enc_finish(e, &itr.itr);
	return e;
}
