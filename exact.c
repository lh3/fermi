#include <assert.h>
#include "rld.h"
#include "kstring.h"
#include "fermi.h"
#include "kvec.h"

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

int64_t fm_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
	uint64_t k = x, *ok;
	ok = alloca(8 * e->asize);
	s->l = 0;
	while (1) {
		int c = rld_rank1a(e, k, ok);
		if (c) {
			k = e->cnt[c] + ok[c] - 1;
			assert(k < e->mcnt[0]);
			kputc(c, s);
		} else return ok[0] - 1;
	}
}

int fm6_extend(const rld_t *e, const fmintv_t *ik, fmintv_t ok[6], int is_back)
{ // FIXME: this can be accelerated a little by using rld_rank1a() when ik.x[2]==1
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

int64_t fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
	fmintv_t ok[6], ik;
	s->l = 0;
	ik.x[0] = ik.x[1] = x; ik.x[2] = 1;
	while (1) {
		int c;
		fm6_extend(e, &ik, ok, 1);
		if (!ok[0].x[2]) {
			for (c = 1; c < 6; ++c)
				if (ok[c].x[2]) break;
			ik = ok[c];
			kputc(c, s);
		} else return ok[0].x[0];
	}
}

/*********************
 *********************/

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

fmintv_t fm6_overlap_intv(const rld_t *e, int len, const uint8_t *seq, int min, int j, int at5, fmintv_v *p)
{ // requirement: seq[j] matches the end of a read
	int c, depth, dir, end;
	fmintv_t ik, ok[6];
	p->n = 0;
	dir = at5? 1 : -1; // at5 is true iff we start from the 5'-end of a read
	end = at5? len : -1;
	c = seq[j];
	fm6_set_intv(e, c, ik);
	for (depth = 1, j += dir; j != end; j += dir, ++depth) {
		c = at5? fm6_comp(seq[j]) : seq[j];
		fm6_extend(e, &ik, ok, !at5);
		if (!ok[c].x[2]) break; // cannot be extended
		if (depth >= min && ok[0].x[2]) {
			//ik.info = j - dir; kv_push(fmintv_t, *p, ik);
			ok[0].info = j - dir; kv_push(fmintv_t, *p, ok[0]);
		}
		ik = ok[c];
	}
	reverse_fmivec(p); // reverse the array such that the smallest interval comes first
	return ik;
}

/****************************
 * SuperMaximal Exact Match *
 ****************************/

int fm6_smem1(const rld_t *e, int len, const uint8_t *q, int x, fmintv_v *mem)
{
	int i, j, c, ret;
	fmintv_t ik, ok[6];
	fmintv_v a[2], *prev, *curr, *swap;

	kv_init(a[0]); kv_init(a[1]);
	prev = &a[0]; curr = &a[1];
	fm6_set_intv(e, q[x], ik);

	ik.info = x + 1;
	for (i = x + 1; i < len; ++i) { // forward search
		c = fm6_comp(q[i]);
		fm6_extend(e, &ik, ok, 0);
		if (ok[c].x[2] != ik.x[2]) { // change of the interval size
			if (ik.x[2] != ok[0].x[2]) kv_push(fmintv_t, *curr, ik);
			if (ok[0].x[2]) { // some sequences come to an end
				ok[0].info = i;
				kv_push(fmintv_t, *curr, ok[0]);
			}
		}
		if (ok[c].x[2] == 0) break; // cannot be extended
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) kv_push(fmintv_t, *curr, ik); // push the last interval if we reach the end
	reverse_fmivec(curr); // s.t. smaller intervals visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;
//	for (i = 0; i < prev->n; ++i) printf("[%lld, %lld, %lld], %lld\n", prev->a[i].x[0], prev->a[i].x[1], prev->a[i].x[2], prev->a[i].info);

	mem->n = 0;
	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			int fl_match; // whether this leads to a full-length read match
			fm6_extend(e, p, ok, 1);
			fl_match = (ok[0].x[2] && p->x[1] < e->mcnt[1]);
			if (ok[c].x[2] == 0 || fl_match || i == -1) { // keep the hit if: full-length match, reaching the beginning or not extended further
				if (curr->n == 0 || fl_match) { // curr->n to make sure there is no longer matches
//					printf("%d, %lld, [%lld,%lld,%lld]\n", i+1, p->info, p->x[0], p->x[1], p->x[2]);
					if (fl_match || mem->n == 0 || i + 1 < (mem->a[mem->n-1].info>>32&FM_MASK30)) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(ok[0].x[2] != 0) << 63 | (uint64_t)(i + 1)<<32; // bit 64 keeps whether the left-end is closed
						kv_push(fmintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			}
			if (ok[c].x[2] && (p->x[1] < e->mcnt[1] || curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = p->info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	reverse_fmivec(mem); // s.t. sorted by the start coordinate

//	for (i = 0; i < mem->n; ++i) printf("[%lld,%lld,%lld] %lld, %lld\n", mem->a[i].x[0], mem->a[i].x[1], mem->a[i].x[2], mem->a[i].info>>32&FM_MASK30, mem->a[i].info&FM_MASK30);
	free(a[0].a); free(a[1].a);
	return ret;
}

int fm6_smem(const rld_t *e, int len, const uint8_t *q, fmintv_v *mem)
{
	int x = 0, i;
	fmintv_v tmp;
	kv_init(tmp);
	mem->n = 0;
	do {
		x = fm6_smem1(e, len, q, x, &tmp);
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
