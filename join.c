#include "fermi.h"
#include "rld.h"
#include "kvec.h"
#include "kstring.h"

int fm6_unambi_nei_for(const rld_t *e, int min, int beg, kstring_t *s, fmintv_v *curr, fmintv_v *prev, uint64_t *bits)
{
	extern fmintv_t fm6_overlap_intv(const rld_t *e, int len, const uint8_t *seq, int min, int j, int at5, fmintv_v *p);
	int i, j, c, old_l = s->l, ret;
	fmintv_t ik, ok[6];
	fmintv_v *swap;

	// backward search for overlapping reads
	ik = fm6_overlap_intv(e, s->l - beg, (uint8_t*)s->s + beg, min, s->l - beg - 1, 0, prev);
	if (prev->n == 0) return -1; // no overlapping reads
	for (j = 0; j < prev->n; ++j) prev->a[j].info += beg;
	ret = prev->a[0].info;
	// test if s[beg..s->l] contained in another read
	fm6_extend(e, &ik, ok, 1);
	if (ik.x[2] != ok[0].x[2]) return -2; // the sequence is left contained
	fm6_extend(e, &ik, ok, 0);
	if (ik.x[2] != ok[0].x[2]) return -2; // the sequence is right contained
	// forward search for forward branching test and for the longest read
	while (prev->n) {
		int c0 = -1;
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			fm6_extend(e, p, ok, 0);
			if (ok[0].x[2]) {
				if ((int)p->info == ret && ok[0].x[2] == p->x[2]) break;
				// then there are contained matches
			}
			if (c0 == -1) {
				for (c = 1; c < 6; ++c) if (ok[c].x[2]) break;
				if (c == 6) continue;
				c0 = c;
			}
			if (ok[c0].x[2] + ok[0].x[2] < p->x[2]) return -3;
			if (ok[c0].x[2] && (curr->n == 0 || ok[c0].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c0].info = p->info;
				kv_push(fmintv_t, *curr, ok[c0]);
			}
		}
		if (j < prev->n) break;
		kputc(fm6_comp(c0), s);
		swap = curr; curr = prev; prev = swap;
	}
	for (i = 0; i < s->l; ++i) putchar("$ACGTN"[(int)s->s[i]]); putchar('\n');

	// forward search for reads overlapping the extension read
	fm6_overlap_intv(e, s->l, (uint8_t*)s->s, min, ret, 1, prev);
	printf("ret=%d, len=%d, prev->n=%d\n", (int)ret, (int)s->l, (int)prev->n);
	// backward search for backward branching test
	for (i = ret - 1; i >= 0 && prev->n; --i) {
		c = s->s[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fm6_extend(e, &prev->a[j], ok, 1);
			if (ok[c].x[2] + ok[0].x[2] != prev->a[j].x[2]) { // branching
				s->l = old_l;
				return -4; // backward branching
			}
			if (ok[c].x[2] && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]))
				kv_push(fmintv_t, *curr, ok[c]);
		}
		swap = curr; curr = prev; prev = swap;
	}
	printf("final: ret=%d, len=%d\n", (int)ret, (int)s->l);
	return ret;
}

void fm6_extend_further1(const rld_t *e, uint64_t x)
{
	extern void seq_reverse(int l, unsigned char *s);
	kstring_t s[2];
	fmintv_v a[2];
	int i;

	kv_init(a[0]); kv_init(a[1]);
	s[0].l = s[0].m = 0; s[0].s = 0;
	s[1].l = s[1].m = 0; s[1].s = 0;
	fm_retrieve(e, x, &s[0]); seq_reverse(s[0].l, (uint8_t*)s[0].s);
	for (i = 0; i < s[0].l; ++i) putchar("$ACGTN"[(int)s[0].s[i]]); putchar('\n');
	printf("*** ret=%d\n", fm6_unambi_nei_for(e, 10, 0, &s[0], &a[0], &a[1], 0));
}

