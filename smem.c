#include "priv.h"
#include "kvec.h"
#include "kstring.h"

int fm6_smem1(const rld_t *e, int len, const uint8_t *q, int x, fmintv_v *mem, int self_match)
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
			if (!self_match && ok[0].x[2]) { // some sequences come to an end
				ok[0].info = i;
				kv_push(fmintv_t, *curr, ok[0]);
			}
		}
		if ((!self_match && ok[c].x[2] == 0) || (self_match && ok[c].x[2] < 2)) break; // cannot be extended
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) kv_push(fmintv_t, *curr, ik); // push the last interval if we reach the end
	fm_reverse_fmivec(curr); // s.t. smaller intervals visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;
//	for (i = 0; i < prev->n; ++i) printf("[%lld, %lld, %lld], %lld\n", prev->a[i].x[0], prev->a[i].x[1], prev->a[i].x[2], prev->a[i].info);

	mem->n = 0;
	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			fmintv_t *p = &prev->a[j];
			int cont, fl_match; // whether this leads to a full-length read match
			fm6_extend(e, p, ok, 1);
			fl_match = (ok[0].x[2] && p->x[1] < e->mcnt[1]);
			cont = self_match? (ok[c].x[2] > 1) : (ok[c].x[2] != 0);
			if (!cont || fl_match || i == -1) { // keep the hit if: full-length match, reaching the beginning or not extended further
				if (curr->n == 0 || fl_match) { // curr->n to make sure there is no longer matches
//					printf("%d, %lld, [%lld,%lld,%lld]\n", i+1, p->info, p->x[0], p->x[1], p->x[2]);
					if (fl_match || mem->n == 0 || i + 1 < (mem->a[mem->n-1].info>>32&FM_MASK30)) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(ok[0].x[2] != 0) << 63 | (uint64_t)(i + 1)<<32; // bit 64 keeps whether the left-end is closed
						kv_push(fmintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			}
			if (cont && (p->x[1] < e->mcnt[1] || curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = p->info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	fm_reverse_fmivec(mem); // s.t. sorted by the start coordinate

//	for (i = 0; i < mem->n; ++i) printf("[%lld,%lld,%lld] %lld, %lld\n", mem->a[i].x[0], mem->a[i].x[1], mem->a[i].x[2], mem->a[i].info>>32&FM_MASK30, mem->a[i].info&FM_MASK30);
	free(a[0].a); free(a[1].a);
	return ret;
}

int fm6_smem(const rld_t *e, int len, const uint8_t *q, fmintv_v *mem, int self_match)
{
	int x = 0, i;
	fmintv_v tmp;
	kv_init(tmp);
	mem->n = 0;
	do {
		x = fm6_smem1(e, len, q, x, &tmp, self_match);
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
