#include <math.h>
#include "priv.h"
#include "kvec.h"
#include "kstring.h"

int fm6_smem1_core(const rld_t *e, int min_occ, int len, const uint8_t *q, int x, fmsmem_v *mem, fmintv_v *prev, fmintv_v *curr)
{ // for more comments, see bwa/bwt.c
	int i, j, c, ret;
	fmintv_t ik, ok[6];
	fmintv_v *swap;

	fm6_set_intv(e, q[x], ik);
	ik.info = x + 1;
	for (i = x + 1; i < len; ++i) { // forward extension
		c = fm6_comp(q[i]);
		fm6_extend(e, &ik, ok, 0);
		if (ok[c].x[2] != ik.x[2]) {
			kv_push(fmintv_t, *curr, ik);
			if (ok[c].x[2] < min_occ) break;
		}
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) kv_push(fmintv_t, *curr, ik);
	kv_reverse(fmintv_t, *curr);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	mem->n = 0;
	for (i = x - 1; i >= -1; --i) {
		c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			fm6_extend(e, p, ok, 1);
			if (c == 0 || ok[c].x[2] < min_occ) {
				if (curr->n == 0) {
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].ik.info>>32) {
						fmsmem_t *q;
						kv_pushp(fmsmem_t, *mem, &q);
						q->ik = *p; q->ik.info |= (uint64_t)(i + 1)<<32;
						memcpy(q->ok[0], ok, 6 * sizeof(fmintv_t));
					}
				}
			} else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
				ok[c].info = p->info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	kv_reverse(fmsmem_t, *mem);
	return ret;
}

int fm6_smem1(const rld_t *e, int len, const uint8_t *q, int x, fmsmem_v *mem, int min_occ)
{
	int ret;
	fmintv_v a[2];
	kv_init(a[0]); kv_init(a[1]);
	ret = fm6_smem1_core(e, min_occ, len, q, x, mem, &a[0], &a[1]);
	free(a[0].a); free(a[1].a);
	return ret;
}

int fm6_smem(const rld_t *e, int len, const uint8_t *q, fmsmem_v *mem, int min_occ)
{
	int x = 0, i;
	fmsmem_v tmp;
	kv_init(tmp);
	mem->n = 0;
	do {
		x = fm6_smem1(e, len, q, x, &tmp, min_occ);
		for (i = 0; i < tmp.n; ++i) {
			kv_push(fmsmem_t, *mem, tmp.a[i]);
		}
	} while (x < len);
	return mem->n;
}

int fm6_write_smem(const rld_t *e, const fmsmem_t *a, kstring_t *s)
{
	s->l = 0;
	kputuw(a->ik.info>>32, s); kputc('\t', s);
	kputuw((uint32_t)a->ik.info, s); kputc('\t', s);
	kputl(a->ik.x[2], s);
	return s->l;
}
