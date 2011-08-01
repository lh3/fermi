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
	k = x + 1;
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

