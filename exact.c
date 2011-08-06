#include "exact.h"

uint64_t fm_backward_search(rld_t *e, const rldidx_t *r, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end)
{
	uint64_t k, l, ok, ol;
	int i, c;
	c = str[len - 1];
	k = e->cnt[c]; l = e->cnt[c + 1] - 1;
	for (i = len - 2; i >= 0; --i) {
		c = str[i];
		rld_rank21(e, r, k - 1, l, c, &ok, &ol);
		k = e->cnt[c] + ok;
		l = e->cnt[c] + ol - 1;
		if (k > l) break;
	}
	if (k > l) return 0;
	*sa_beg = k; *sa_end = l;
	return l - k + 1;
}

void fm_retrieve(rld_t *e, const rldidx_t *r, uint64_t x, kstring_t *s)
{
	uint64_t k, *ok, *ol;
	ok = alloca(8 * e->asize);
	ol = alloca(8 * e->asize);
	s->l = 0;
	k = x;
	while (1) {
		int c;
		rld_rank2a(e, r, k - 1, k, ok, ol);
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

int fm6_extend(rld_t *e, const rldidx_t *r, uint64_t ik[3], uint64_t ok[24], int is_back)
{
	uint64_t tk[6], tl[6];
	int i;
	rld_rank2a(e, r, ik[!is_back] - 1, ik[!is_back] - 1 + ik[2], tk, tl);
	for (i = 0; i < 6; ++i) {
		ok[i<<2|(!is_back)] = e->cnt[i] + tk[i];
		ok[i<<2|2] = (tl[i] -= tk[i]);
	}
	ok[is_back+0]  = ik[is_back+0];
	ok[is_back+4]  = ok[is_back+0]  + tl[0];
	ok[is_back+8]  = ok[is_back+4]  + tl[4];
	ok[is_back+12] = ok[is_back+8]  + tl[3];
	ok[is_back+16] = ok[is_back+12] + tl[2];
	ok[is_back+20] = ok[is_back+16] + tl[1];
	return 0;
}

int fm6_extend1(rld_t *e, const rldidx_t *r, uint64_t ik[3], int c, uint64_t ok[3], int is_back)
{
	uint64_t o[24];
	if (!is_back) c = (c >= 1 && c <= 4)? 5 - c : c;
	fm6_extend(e, r, ik, o, is_back);
	ok[0] = o[c<<2|0]; ok[1] = o[c<<2|1]; ok[2] = o[c<<2|2];
	return 0;
}

void fm6_retrieve(rld_t *e, const rldidx_t *r, uint64_t x, kstring_t *s)
{
	uint64_t ok[24], ik[3];
	s->l = 0;
	ik[0] = ik[1] = x; ik[2] = 1;
	while (1) {
		int c;
		fm6_extend(e, r, ik, ok, 0);
		for (c = 0; c < 6; ++c)
			if (ok[c<<2|2] == 1) break;
		/*{
			int i, j;
			for (i = 0; i < 6; ++i) {
				for (j = 0; j < 3; ++j)
					printf("%2lld, ", ok[i<<2|j]);
				printf("| ");
			}
			putchar('\n');
		}*/
		if (c) {
			ik[0] = ok[c<<2], ik[1] = ok[c<<2|1], ik[2] = ok[c<<2|2];
			kputc(c, s);
		} else break;
	}
}
