#include <stdio.h>
#include "fermi.h"
#include "kstring.h"

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
