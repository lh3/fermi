#include "fermi.h"
#include "rld.h"
#include "kstring.h"

void msg_write_node(const fmnode_t *p, long id, kstring_t *out)
{
	int j, k;
	kputc('@', out); kputl(id, out);
	for (j = 0; j < 2; ++j) {
		kputc('\t', out);
		kputl(p->k[j], out); kputc('>', out);
		for (k = 0; k < p->nei[j].n; ++k) {
			if (k) kputc(',', out);
			kputl(p->nei[j].a[k].x, out); kputc(':', out);
			kputw((int32_t)p->nei[j].a[k].y, out); kputc(':', out);
			kputw((int32_t)(p->nei[j].a[k].y>>32), out);
		}
		if (p->nei[j].n == 0) kputc('.', out);
	}
	kputc('\n', out);
	ks_resize(out, out->l + 2 * p->l + 5);
	for (j = 0; j < p->l; ++j) out->s[out->l++] = "ACGT"[(int)p->seq[j] - 1];
	out->s[out->l] = 0;
	kputsn("\n+\n", 3, out);
	kputsn(p->cov, p->l, out);
	kputc('\n', out);
}

void msg_print(const fmnode_v *nodes)
{
	kstring_t out;
	size_t i;
	out.l = out.m = 0; out.s = 0;
	for (i = 0; i < nodes->n; ++i) {
		out.l = 0;
		msg_write_node(&nodes->a[i], i, &out);
		fputs(out.s, stdout);
	}
	free(out.s);
}
