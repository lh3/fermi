#include <stdio.h>
#include "fermi.h"
#include "kstring.h"
#include "kvec.h"

typedef kvec_t(fmgelem_t*) fmgelem_pv;

#include "khash.h"
KHASH_MAP_INIT_INT64(g, fmgelem_pv)

static void add(khash_t(g) *h, uint64_t x, fmgelem_t *a)
{
	khint_t k;
	int ret;
	fmgelem_pv *p;
	if (x == 0) return;
	k = kh_put(g, h, x, &ret);
	p = &kh_val(h, k);
	if (ret) kv_init(*p);
	kv_push(fmgelem_t*, *p, a);
}

void bog_output(const fmgelem_v *elems)
{
	size_t i;
	kstring_t out;
	int j;
	out.l = out.m = 0; out.s = 0;
	for (i = 0; i < elems->n; ++i) {
		fmgelem_t *p = &elems->a[i];
		printf("@%ld %c%c%ld %c%c%ld\n", i, p->type[0], p->dir[0], (long)p->k[0], p->type[1], p->dir[1], (long)p->k[1]);
		out.l = 0;
		ks_resize(&out, p->l + 1);
		for (j = 0; j < p->l; ++j) out.s[out.l++] = "ACGT"[(int)p->seq[j]];
		out.s[out.l] = 0;
		puts(out.s);
		puts("+");
		memcpy(out.s, p->cov, p->l);
		puts(out.s);
	}
}

void bog_clean(fmgelem_v *elems)
{
	size_t i;
	khash_t(g) *h;
	khint_t k;
	h = kh_init(g);
	for (i = 0; i < elems->n; ++i) {
		add(h, elems->a[i].k[0], &elems->a[i]);
		if (elems->a[i].k[0] != elems->a[i].k[1])
			add(h, elems->a[i].k[1], &elems->a[i]);
	}
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			if (kh_val(h, k).n > 1) fprintf(stderr, "x\t%ld\t%ld\n", (long)kh_key(h, k), kh_val(h, k).n);
		}
	}
}
