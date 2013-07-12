#include "priv.h"
#include "kvec.h"

void fmc_opt_init(fmec2opt_t *opt)
{
	opt->min_l = 17;
	opt->min_occ = 6;
	opt->min_occ_patch = 3;
	opt->len_factor = 3;

	opt->qual_plus = 10;
	opt->max_pen = 60;
	opt->max_d = 7;
	opt->diff_factor = 13;
	opt->ratio_factor = 10.;
}

void fmc_aux_destroy(fmec2aux_t *a)
{
	free(a->tmp[0].a); free(a->tmp[1].a); free(a->mem1.a); free(a->mem.a);
	free(a->seq); free(a->s); free(a->matrix);
	free(a->f.a); free(a->b.a);
	free(a);
}

static inline void ec_cns_gen(const fmec2opt_t *opt, const fmintv_t k[6], uint8_t q[4])
{
	int c;
	int64_t sum = 0;
	for (c = 1; c <= 4; ++c) sum += k[c].x[2];
	for (c = 1; c <= 4; ++c) {
		int64_t d0, d1;
		d0 = k[c].x[2]; d1 = sum - d0;
		if (d0 < d1) {
			int p, r;
			p = ((d1 < opt->max_d? d1 : opt->max_d) - d0) * opt->diff_factor;
			p = p > 1? p : 1;
			p = p < opt->max_pen? p : opt->max_pen;
			r = d0? (int)(opt->ratio_factor * d1 / d0 + .499) : opt->max_pen;
			q[c-1] = p < r? p : r;
		} else q[c-1] = 0;
	}
}

static long ec_ecinfo(const fmec2opt_t *opt, const fmintv_t k[6], int l_seq, const uint8_t *seq, const fmec2seq_t *s, int pos, int is_comp)
{
	int c;
	uint8_t q[4];
	ec_cns_gen(opt, k, q);
	if (pos >= 0 && pos < l_seq) { // FIXME: more carefully deal with triallelic sites!
		int min, min2, min_c = -1, min2_c = -1, b_read, b_cns, q_cns, penalty;
		b_read = is_comp? fm6_comp(seq[pos]) : seq[pos];
		for (c = 0; c < 4; ++c)
			if (c != b_read - 1)
				q[c] += s[pos].oq;
		for (c = 0, min = min2 = 256; c < 4; ++c) {
			if (min > q[c]) min2 = min, min2_c = min_c, min = q[c], min_c = c;
			else if (min2 > q[c]) min2 = q[c], min2_c = c;
		}
		b_cns = is_comp? 3 - min_c : min_c; // top consensus base
		q_cns = min2 - min; // consensus quality
		q_cns = q_cns < opt->max_pen? q_cns : opt->max_pen;
		penalty = b_read < 5? q[b_read-1] - min : 0;
		penalty = penalty < opt->max_pen? penalty : opt->max_pen;
		return b_cns | q_cns << 4 | (long)penalty << 16;
	} else return 0;
}

static void ec_patch(const fmec2opt_t *opt, const rld_t *e, fmintv_t k[6], int beg, int end, int l_seq, const uint8_t *seq, fmec2seq_t *s)
{
	int i, is_back = (beg > end), step = is_back? -1 : 1;
	if (fm_verbose >= 5) fprintf(stderr, "===> Patching region [%d, %d) <===\n", beg, end);
	for (i = beg; i != end; i += step) {
		int c;
		s[i].cb[is_back] = (k[0].info & 0xf) + 1;
		s[i].cq[is_back] = k[0].info >> 4 & 0xff;
		s[i].cf |= 2<<is_back;
		if (fm_verbose >= 6 || (fm_verbose >= 5 && seq[i] != s[i].cb[is_back])) {
			int c;
			fprintf(stderr, "pos=%d, %c%d => %c%d, depth=(%ld", i, "$ACGTN"[seq[i]], s[i].oq, "$ACGTN"[s[i].cb[is_back]], s[i].cq[is_back], (long)k[0].x[2]);
			if (is_back) for (c = 1; c < 5; ++c) fprintf(stderr, ",%ld", (long)k[c].x[2]);
			else for (c = 4; c > 0; --c) fprintf(stderr, ",%ld", (long)k[c].x[2]);
			fprintf(stderr, ",%ld)\n", (long)k[5].x[2]);
		}
		c = is_back? s[i].cb[is_back] : fm6_comp(s[i].cb[is_back]);
		if (i == end + step || k[c].x[2] < opt->min_occ_patch) break;
		fm6_extend(e, &k[c], k, is_back);
		k[0].info = ec_ecinfo(opt, k, l_seq, seq, s, i + step, !is_back);
	}
}

int fmc_ec_core(const fmec2opt_t *opt, const rld_t *e, fmec2aux_t *aux, int l_seq, char *seq, char *qual)
{
	int i, n_row, x = 0, n_corr, l_cov;
	// allocate enough memory
	if (l_seq > aux->max_len) {
		aux->max_len = l_seq;
		kroundup32(aux->max_len);
		aux->s = realloc(aux->s, aux->max_len * sizeof(fmec2seq_t));
		aux->seq = realloc(aux->seq, aux->max_len);
	}
	// change encoding
	memset(aux->s, 0, l_seq * sizeof(fmec2seq_t));
	for (i = 0; i < l_seq; ++i) {
		aux->seq[i] = seq[i] >= 0? seq_nt6_table[(int)seq[i]] : 5;
		aux->s[i].oq = qual[i] - 33;
	}
	// fill the aux->mem vector
	aux->mem.n = 0;
	do {
		x = fm6_smem1_core(e, opt->min_occ, l_seq, aux->seq, x, &aux->mem1, &aux->tmp[0], &aux->tmp[1]);
		for (i = 0; i < aux->mem1.n; ++i) {
			fmsmem_t *p = &aux->mem1.a[i];
			if ((uint32_t)p->ik.info - (p->ik.info>>32) >= opt->min_l) {
				fm6_extend(e, &p->ik, p->ok[1], 0);
				kv_push(fmsmem_t, aux->mem, *p);
			}
		}
	} while (x < l_seq);
	// make sure the memory allocated to matrix is large enough
	n_row = aux->mem.n + 2;
	x = n_row * n_row;
	if (x > aux->max_matrix) {
		aux->max_matrix = x;
		kroundup32(aux->max_matrix);
		aux->matrix = calloc(aux->max_matrix, sizeof(int));
	}
	// compute the consensus at the ends of segments
	for (i = 0; i < aux->mem.n; ++i) {
		fmsmem_t *p = &aux->mem.a[i];
		p->ok[0][0].info = ec_ecinfo(opt, p->ok[0], l_seq, aux->seq, aux->s, (int)(p->ik.info>>32) - 1, 0);
		p->ok[1][0].info = ec_ecinfo(opt, p->ok[1], l_seq, aux->seq, aux->s, (int32_t)p->ik.info, 1);
	}
	// fill the scoring matrix
	memset(aux->matrix, 0, n_row * n_row * sizeof(int));
	for (i = 0; i < aux->mem.n; ++i) { // first row
		fmsmem_t *p = &aux->mem.a[i];
		aux->matrix[i + 1] = opt->len_factor * ((int32_t)p->ik.info - (p->ik.info>>32));
	}
	for (i = 0; i < aux->mem.n; ++i) { // second to (n_row-2)-th row
		fmsmem_t *pi = &aux->mem.a[i];
		int j, *mat = &aux->matrix[(i + 1) * n_row];
		for (j = i + 1; j < aux->mem.n; ++j) {
			fmsmem_t *pj = &aux->mem.a[j];
			if (pj->ik.info>>32 <= (int32_t)pi->ik.info) { // has overlap
				int l = (int32_t)pj->ik.info - (int32_t)pi->ik.info - 2;
				mat[j+1] = l * opt->len_factor - (int)(pi->ok[1][0].info>>16) - (int)(pj->ok[0][0].info>>16);
			} else { // no overlap
				mat[j+1] = ((int32_t)pj->ik.info - (pj->ik.info>>32)) * opt->len_factor;
			}
		}
	}
	// dynamic programming
	kv_resize(int32_t, aux->f, n_row);
	kv_resize(int32_t, aux->b, n_row);
	aux->f.n = aux->b.n = n_row;
	aux->f.a[0] = aux->b.a[0] = 0;
	for (i = 1; i < n_row; ++i) {
		int j, max, max_j, *mat = &aux->matrix[(i-1) * n_row + i];
		for (j = i - 1, max_j = -1, max = -0x3fffffff; j >= 0; --j, mat -= n_row) {
			int x = aux->f.a[j] + (*mat); // (*mat) = aux->matrix[j * n_row + i]
			if (max < x) max = x, max_j = j;
		}
		aux->f.a[i] = max; aux->b.a[i] = max_j;
	}
	// backtrack
	aux->f.n = 0;
	i = aux->b.a[n_row-1];
	while (i) {
		kv_push(int, aux->f, i - 1);
		i = aux->b.a[i];
	}
	kv_reverse(int, aux->f);
	// debugging information
	if (fm_verbose >= 5) {
		fprintf(stderr, "===> %ld segments <===\n", aux->mem.n);
		for (i = 0; i < aux->mem.n; ++i) {
			fmsmem_t *p = &aux->mem.a[i];
			int j, *mat = &aux->matrix[(i+1) * n_row], beg = p->ik.info>>32, end = (uint32_t)p->ik.info;
			fprintf(stderr, "%4d%4d%4ld  %c%c |", beg, end, (long)p->ik.x[2], "$ACGTN"[beg?aux->seq[beg-1]:0], "$ACGTN"[end<l_seq?aux->seq[end]:0]);
			for (j = 1; j < 5; ++j)
				fprintf(stderr, " %c:%-2ld", "$ACGTN"[j], (long)p->ok[0][j].x[2]);
			fprintf(stderr, " |");
			for (j = 1; j < 5; ++j)
				fprintf(stderr, " %c:%-2ld", "$ACGTN"[j], (long)p->ok[1][j].x[2]);
			fprintf(stderr, " |%3ld%3ld |", (long)p->ok[0][0].info>>16, (long)p->ok[1][0].info>>16);
			for (j = i + 1; j < aux->mem.n; ++j)
				fprintf(stderr, "%5d", mat[j+1]);
			fputc('\n', stderr);
		}
		fprintf(stderr, "===> %ld segments selected <===\n", aux->f.n);
		for (i = 0; i < aux->f.n; ++i) {
			fmsmem_t *p = &aux->mem.a[aux->f.a[i]];
			fprintf(stderr, "%d\t%d\n", (int)(p->ik.info>>32), (uint32_t)p->ik.info);
		}
	}
	// propose corrections
	aux->mem.n = 0;
	for (i = 0; i < aux->f.n; ++i)
		kv_push(fmsmem_t, aux->mem, aux->mem.a[aux->f.a[i]]);
	for (i = 0; i < aux->mem.n; ++i) {
		fmsmem_t *p = &aux->mem.a[i];
		int j, prev, next, beg = p->ik.info>>32, end = (uint32_t)p->ik.info;
		for (j = beg; j < end; ++j) aux->s[j].cf |= 1;
		prev = i? (uint32_t)(p-1)->ik.info : 0;
		next = i < aux->mem.n-1? (p+1)->ik.info>>32 : l_seq;
		if (prev < beg) ec_patch(opt, e, p->ok[0], beg-1, prev-1, l_seq, aux->seq, aux->s);
		if (next > end) ec_patch(opt, e, p->ok[1], end, next, l_seq, aux->seq, aux->s);
	}
	// correct errors
	for (i = 0, n_corr = l_cov = 0; i < l_seq; ++i) {
		fmec2seq_t *si = &aux->s[i];
		int c = aux->seq[i], q = si->oq>>1<<1;
		if ((si->cf>>1&3) == 3) { // inspected from both strands
			if (si->cb[0] != si->cb[1]) {
				if (si->cq[0] > si->cq[1]) c = si->cb[0], q = si->cq[0] - si->cq[1];
				else c = si->cb[1], q = si->cq[1] - si->cq[0];
			} else c = si->cb[0], q = si->cq[0] > si->cq[1]? si->cq[0] : si->cq[1];
			q = q>>1<<1 | 1;
			++l_cov;
		} else if (si->cf>>1&1) { // inspected from one strand
			c = si->cb[0], q = si->cq[0]>>1<<1 | 1;
			++l_cov;
		} else if (si->cf>>1&2) {
			c = si->cb[1], q = si->cq[1]>>1<<1 | 1;
			++l_cov;
		} else if (si->cf&1) { // covered by a SMEM
			q = q + opt->qual_plus < opt->max_pen? q + opt->qual_plus : opt->max_pen;
			q = q>>1<<1 | 1;
			++l_cov;
		}
		n_corr += (aux->seq[i] != c);
		seq[i] = aux->seq[i] == c? "$ACGTN"[c] : "$acgtn"[c];
		qual[i] = q + 33;
	}
	return l_cov;
}
