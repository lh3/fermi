#include <assert.h>
#include "rld.h"
#include "kstring.h"
#include "fermi.h"
#include "kvec.h"

#define fm6_set_intv(e, c, ik) ((ik).x[0] = (e)->cnt[(c)], (ik).x[2] = (e)->cnt[c+1] - (e)->cnt[c], (ik).x[1] = (e)->cnt[fm6_comp(c)])

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
	uint64_t k, *ok;
	ok = alloca(8 * e->asize);
	s->l = 0;
	k = x;
	while (1) {
		int c = rld_rank1a(e, k, ok);
		if (c) {
			k = e->cnt[c] + ok[c] - 1;
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

/*********************
 *********************/

typedef kvec_t(fmintv_t) fmintv_v;
typedef kvec_t(uint64_t) vec64_t;

static inline void reverse_fmivec(fmintv_v *p)
{
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			fmintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}

static inline void overlap_intv(const rld_t *e, int len, const uint8_t *seq, int min, int j, int at5, fmintv_v *p)
{ // requirement: seq[j] matches the end of a read
	int c, depth, dir, end;
	fmintv_t ik, ok[6];
	p->n = 0;
	dir = at5? 1 : -1; // at5 is true iff we start from the 5'-end of a read
	end = at5? len - 1 : 0;
	c = seq[j];
	fm6_set_intv(e, c, ik);
	for (depth = 1, j += dir; j != end; j += dir, ++depth) {
		c = at5? fm6_comp(seq[j]) : seq[j];
		fm6_extend(e, &ik, ok, !at5);
		if (!ok[c].x[2]) break; // cannot be extended
		if (depth >= min && ok[0].x[2] && (p->n == 0 || ik.x[2] != p->a[p->n-1].x[2]))
			kv_push(fmintv_t, *p, ik);
		ik = ok[c];
	}
	reverse_fmivec(p); // reverse the array such that smaller intervals come first
}

int fm6_search_overlap(const rld_t *e, int min, int len, const uint8_t *seq, int is_back)
{
	int i, j, c, last_sentinel, last_beg, shift;
	fmintv_t ik, ok[6];
	fmintv_v a[2], *curr, *prev, *tmp;

	// initialize the vectors
	kv_init(a[0]); kv_init(a[1]);
	prev = &a[0]; curr = &a[1];
	fm6_set_intv(e, seq[is_back? len-1 : 0], ik);
	kv_push(fmintv_t, *prev, ik);

	// core loop
	if (!is_back) last_sentinel = last_beg = 0,       shift = -1, i = 1;
	else          last_sentinel = last_beg = len - 1, shift =  1, i = len - 2;
	for (;;) {
		c = is_back? seq[i] : fm6_comp(seq[i]);
		curr->n = 0; // clear the curr list
		for (j = 0; j < prev->n; ++j) { // traverse the prev list
			fm6_extend(e, &prev->a[j], ok, is_back);
			if (ok[c].x[2] && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]))
				kv_push(fmintv_t, *curr, ok[c]);
			if (ok[0].x[2]) last_sentinel = i + shift;
		}
		if (curr->n == 0) { // backward search
			if (!is_back) {
				if (last_sentinel < last_beg) break;
			} else {
				if (last_sentinel > last_beg) break;
			}
			overlap_intv(e, len, seq, min, last_sentinel, is_back, curr);
			if (curr->n == 0) break;
			i = last_sentinel; // i will be increased by 1 in the next round of the loop
			last_beg = i - shift;
		}
		tmp = curr; curr = prev; prev = tmp;
		if (!is_back) {
			if (++i >= len) break;
		} else {
			if (--i < 0) break;
		}
	}
	kv_destroy(a[0]); kv_destroy(a[1]);
	return is_back? len - 1 - i : i;
}
