#include "rld.h"
#include "kstring.h"
#include "fermi.h"

uint64_t fm_backward_search(const rld_t *e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end)
{
	uint64_t k, l, ok, ol;
	int i, c;
	c = str[len - 1];
	k = e->cnt[c]; l = e->cnt[c + 1] - 1;
	for (i = len - 2; i >= 0; --i) {
		c = str[i];
		rld_rank21(e, k - 1, l, c, &ok, &ol);
		k = e->cnt[c] + ok;
		l = e->cnt[c] + ol - 1;
		if (k > l) break;
	}
	if (k > l) return 0;
	*sa_beg = k; *sa_end = l;
	return l - k + 1;
}

void fm_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
	uint64_t k, *ok, *ol;
	ok = alloca(8 * e->asize);
	ol = alloca(8 * e->asize);
	s->l = 0;
	k = x;
	while (1) {
		int c;
		rld_rank2a(e, k - 1, k, ok, ol);
		for (c = 0; c < e->asize; ++c) { // FIXME: to simplify
			ok[c] += e->cnt[c];
			ol[c] += e->cnt[c] - 1;
			if (ol[c] == ok[c]) break;
		}
		if (c) {
			k = ok[c];
			kputc(c, s);
		} else break;
	}
}

int fm6_extend(const rld_t *e, const fmintv_t *ik, fmintv_t ok[6], int is_back)
{
	uint64_t tk[6], tl[6];
	int i;
	rld_rank2a(e, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i < 6; ++i) {
		ok[i].x[!is_back] = e->cnt[i] + tk[i];
		ok[i].x[2] = (tl[i] -= tk[i]);
	}
	ok[0].x[is_back] = ik->x[is_back];
	ok[4].x[is_back] = ok[0].x[is_back] + tl[0];
	ok[3].x[is_back] = ok[4].x[is_back] + tl[4];
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[5].x[is_back] = ok[1].x[is_back] + tl[1];
	return 0;
}

void fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
	fmintv_t ok[6], ik;
	s->l = 0;
	ik.x[0] = ik.x[1] = x; ik.x[2] = 1;
	while (1) {
		int c;
		fm6_extend(e, &ik, ok, 0);
		for (c = 0; c < 6; ++c)
			if (ok[c].x[2] == 1) break;
		if (c) {
			ik.x[0] = ok[c].x[0], ik.x[1] = ok[c].x[1], ik.x[2] = ok[c].x[2];
			kputc(c, s);
		} else break;
	}
}

/*
#include "khash.h"

#define fmintv_eq(a, b) ((a).x[0] == (b).x[0] && (a).x[2] == (b).x[2])
#define fmintv_hash(a) ((a).x[0] ^ (a).x[2]<<11)
KHASH_INIT2(fm,, fmintv_t, char, 0, fmintv_hash, fmintv_eq)
*/

#include "kvec.h"
#include "ksort.h"
//#define intvcmp(a, b) ((a).x[2] < (b).x[2] || ((a).x[2] == (b).x[2] && (a).x[1] < (b).x[1]))
#define intvcmp(a, b) ((a).x[2] < (b).x[2])
KSORT_INIT(fm, fmintv_t, intvcmp)

#define fm6_set_intv(e, c, ik) ((ik).x[0] = (e)->cnt[(c)], (ik).x[2] = (e)->cnt[c+1] - (e)->cnt[c], (ik).x[1] = (e)->cnt[fm6_comp(c)])

#include <stdio.h>

int fm6_search_forward_overlap(const rld_t *e, int min, int len, const uint8_t *seq)
{
	int i, j, k, c, last_sentinel = 0, last_beg = 0;
	fmintv_t ik, ok[6];
	kvec_t(fmintv_t) a[2], *curr, *prev, *tmp;

	// initialize the vectors
	kv_init(a[0]); kv_init(a[1]);
	prev = &a[0]; curr = &a[1];
	fm6_set_intv(e, seq[0], ik);
	kv_push(fmintv_t, *prev, ik);

	// core loop
	for (i = 1; i < len; ++i) {
		c = fm6_comp(seq[i]);
		curr->n = 0; // clear the curr list
		for (j = 0; j < prev->n; ++j) { // traverse the prev list
			fm6_extend(e, &prev->a[j], ok, 0);
			if (ok[c].x[2]) kv_push(fmintv_t, *curr, ok[c]);
			if (ok[0].x[2]) last_sentinel = i - 1;
		}
		if (curr->n == 0) {
			int depth;
			if (last_sentinel < last_beg) break;
			c = seq[last_sentinel];
			fm6_set_intv(e, c, ik);
			for (j = last_sentinel - 1, depth = 1; j > 0; --j, ++depth) {
				c = seq[j];
				fm6_extend(e, &ik, ok, 1);
				if (!ok[c].x[2]) break;
				if (depth >= min && ok[0].x[2]) kv_push(fmintv_t, *curr, ik);
				ik = ok[c];
			}
			if (curr->n == 0) break;
			i = last_sentinel; // i will be increased by 1 in the next round of the loop
			last_beg = i + 1;
		}
		if (curr->n > 1) { // then de-redundancy
			ks_introsort(fm, curr->n, curr->a); // FIXME: reconsider if sorting is necessary; maybe not
			for (j = k = 1; j < curr->n; ++j) {
				if (curr->a[k].x[2] != curr->a[j].x[2]) {
					if (k != j) curr->a[k++] = curr->a[j];
					else ++k;
				}
			}
			curr->n = k;
		}
		tmp = curr; curr = prev; prev = tmp;
	}
	kv_destroy(a[0]); kv_destroy(a[1]);
	return i;
}
