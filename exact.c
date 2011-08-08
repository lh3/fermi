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

int fm6_extend(const rld_t *e, fmintv_t ik, fmintv_t ok[6], int is_back)
{
	uint64_t tk[6], tl[6];
	int i;
	rld_rank2a(e, ik.x[!is_back] - 1, ik.x[!is_back] - 1 + ik.x[2], tk, tl);
	for (i = 0; i < 6; ++i) {
		ok[i].x[!is_back] = e->cnt[i] + tk[i];
		ok[i].x[2] = (tl[i] -= tk[i]);
	}
	ok[0].x[is_back] = ik.x[is_back];
	ok[1].x[is_back] = ok[0].x[is_back] + tl[0];
	ok[2].x[is_back] = ok[1].x[is_back] + tl[4];
	ok[3].x[is_back] = ok[2].x[is_back] + tl[3];
	ok[4].x[is_back] = ok[3].x[is_back] + tl[2];
	ok[5].x[is_back] = ok[4].x[is_back] + tl[1];
	return 0;
}

void fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
	fmintv_t ok[6], ik;
	s->l = 0;
	ik.x[0] = ik.x[1] = x; ik.x[2] = 1;
	while (1) {
		int c;
		fm6_extend(e, ik, ok, 0);
		for (c = 0; c < 6; ++c)
			if (ok[c].x[2] == 1) break;
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
			ik.x[0] = ok[c].x[0], ik.x[1] = ok[c].x[1], ik.x[2] = ok[c].x[2];
			kputc(c, s);
		} else break;
	}
}
