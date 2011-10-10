#include "fermi.h"
#include "rld.h"
#include "kstring.h"

void msg_print(const fmnode_v *nodes)
{
	kstring_t out;
	size_t i;
	int j, k;
	out.l = out.m = 0; out.s = 0;
	for (i = 0; i < nodes->n; ++i) {
		fmnode_t *p = &nodes->a[i];
		out.l = 0;
		kputc('@', &out); kputl(i, &out);
		for (j = 0; j < 2; ++j) {
			kputc(' ', &out);
			kputl(p->k[j], &out); kputc('>', &out);
			for (k = 0; k < p->nei[j].n; ++k) {
				if (k) kputc(',', &out);
				kputl(p->nei[j].a[k], &out);
			}
			if (p->nei[j].n == 0) kputc('.', &out);
		}
		kputc('\n', &out);
		ks_resize(&out, out.l + p->l + 3);
		for (j = 0; j < p->l; ++j) out.s[out.l++] = "ACGT"[(int)p->seq[j] - 1];
		kputs("\n+", &out);
		out.s[out.l] = 0;
		puts(out.s);
		puts(p->cov);
	}
}
