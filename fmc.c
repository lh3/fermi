#include <zlib.h>
#include <string.h>
#include <pthread.h>
#include "priv.h"
#include "kstring.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

typedef struct {
	int n_threads;
	int min_l, min_occ;
	int for_qthres, rev_qthres;
	int multi_thres;
	int chunk_size;
	int diff_factor;
	double ratio_factor;
} fmcopt_t;

typedef struct {
	fmintv_t k[6];
} fmintv6_t;

typedef struct { size_t n, m; fmintv6_t *a; } fmintv6_v;

typedef struct {
	uint32_t b0:3, b1:3, ec:2, q0:8, q1:8;
} ecrst_t;

void fmc_opt_init(fmcopt_t *opt)
{
	memset(opt, 0, sizeof(fmcopt_t));
	opt->n_threads = 1;
	opt->min_l = 23;
	opt->min_occ = 4;
	opt->for_qthres = 30;
	opt->rev_qthres = 20;
	opt->multi_thres = 10;
	opt->diff_factor = 13;
	opt->chunk_size = 1<<26;
	opt->ratio_factor = .8;
}

static inline ecrst_t ec_recommend(const fmcopt_t *opt, const fmintv_t k[6], int base, int qual)
{
	int c, q, b0, b1;
	int64_t sum, max0, max1;
	ecrst_t r;

	r.b0 = base; r.ec = 0; r.b1 = 0; r.q0 = r.q1 = 0;
	b0 = b1 = 0; max0 = max1 = sum = 0;
	for (c = 1; c <= 4; ++c) {
		if (k[c].x[2] > max0) max1 = max0, b1 = b0, max0 = k[c].x[2], b0 = c;
		else if (k[c].x[2] > max1) max1 = k[c].x[2], b1 = c;
		sum += k[c].x[2];
	}
	max1 += max0;
	if (max0 == sum) { // no other types of bases
		q = max0 * opt->diff_factor;
		q = q < 60? q : 60;
		if (base == b0) r.q0 = qual > q? qual : q;
		else if (q <= qual) r.q0 = qual - q;
		else r.b0 = b0, r.q0 = q - qual, r.ec = 1;
	} else { // in the following, max0 < sum
		int tmp, q1;
		// q0
		q = (max0 - (sum - max0)) * opt->diff_factor;
		tmp = (7 - (sum - max0)) * opt->diff_factor;
		q = q < tmp? q : tmp;
		tmp = (int)(opt->ratio_factor * max0 / (sum - max0) + .499);
		q = q < tmp? q : tmp;
		q = q > 1? q : 1;
		if (base != b0) {
			if (base != b1) {
				q1 = (max1 - (sum - max1)) * opt->diff_factor;
				tmp = (7 - (sum - max1)) * opt->diff_factor;
				q1 = q1 < tmp? q1 : tmp;
				tmp = sum != max1? (int)(opt->ratio_factor * max1 / (sum - max1) + .499) : 60;
				q1 = q1 < tmp? q1 : tmp;
				q1 = q1 > 1? q1 : 1;
				if (q1 > qual) {
					r.b0 = b0; r.q0 = q1 - qual;
					r.b1 = b1; r.q1 = q;
					r.ec = 2;
				} else r.q0 = qual - q1;
			} else {
				if (q <= qual) r.q0 = qual - q;
				else r.b0 = b0, r.q0 = q - qual, r.ec = 1;
			}
		} else r.q0 = qual > q? qual : q;
	}
	return r;
}

int fm6_ec2_core(const fmcopt_t *opt, const rld_t *e, int l_seq, uint8_t *seq, uint8_t *qual, int x, fmintv_v *prev, fmintv_v *curr, fmintv6_v *tmp)
{ // this function is similar to fm6_smem1_core()
	int i, j, c, ret;
	fmintv_t ik;
	fmintv_v *swap;
	fmintv6_t ok;
	
	fm6_set_intv(e, seq[x], ik);
	ik.info = x + 1;

	for (i = x + 1; i < l_seq; ++i) { // forward search
		c = fm6_comp(seq[i]);
		fm6_extend(e, &ik, ok.k, 0); // forward extension
		if (ok.k[c].x[2] != ik.x[2]) { // change of interval size
			ecrst_t r;
			kv_push(fmintv_t, *curr, ik);
			if (i - x >= opt->min_l) {
				r = ec_recommend(opt, ok.k, c, qual[i]);
				if (r.ec == 1 && r.q0 > opt->for_qthres)
					c = r.b0, seq[i] = fm6_comp(c);
			}
			if (ok.k[c].x[2] < opt->min_occ) break;
			// TODO: should I consider ambiguous bases more carefully?
		}
		ik = ok.k[c]; ik.info = i + 1;
	}

	if (i == l_seq) kv_push(fmintv_t, *curr, ik);
	fm_reverse_fmivec(curr);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= 0; --i) { // backward search
		int cc, cq;
		cc = c = seq[i];
		kv_resize(fmintv6_t, *tmp, prev->n);
		tmp->n = prev->n;
		for (j = 0; j < prev->n; ++j) // collect all the intervals at the current position
			fm6_extend(e, &prev->a[j], tmp->a[j].k, 1);
		cc = 0; cq = 0;
		for (j = 0; j < prev->n; ++j) { // check if we need to make a correction
			if (prev->a[j].info - i >= opt->min_l) {
				ecrst_t r;
				r = ec_recommend(opt, tmp->a[j].k, c, qual[i]);
				if (r.ec == 1 && cq < r.q0) cc = r.b0, cq = r.q0;
			}
		}
		if (cq > opt->rev_qthres) seq[i] = c = cc;
		for (j = 0; j < prev->n; ++j) { // update $curr
			fmintv_t *ok = tmp->a[j].k;
			if (ok[c].x[2] >= opt->min_occ && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = prev->a[j].info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	return ret;
}

static int fm6_ec2_corr1(const fmcopt_t *opt, const rld_t *e, int l_seq, uint8_t *seq, uint8_t *qual, fmintv_v *prev, fmintv_v *curr, fmintv6_v *tmp)
{
	int x = 0;
	while ((x = fm6_ec2_core(opt, e, l_seq, seq, qual, x, prev, curr, tmp)) < l_seq);
	seq_revcomp6(l_seq, seq);
	seq_reverse(l_seq, qual);
	x = 0;
	while ((x = fm6_ec2_core(opt, e, l_seq, seq, qual, x, prev, curr, tmp)) < l_seq);
	seq_revcomp6(l_seq, seq);
	seq_reverse(l_seq, qual);
	return 0;
}

typedef struct {
	int n, step, start;
	const fmcopt_t *opt;
	const rld_t *e;
	const int *len;
	uint8_t **seq;
	uint8_t **qual;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int i, j;
	fmintv_v prev, curr;
	fmintv6_v tmp;

	memset(&prev, 0, sizeof(fmintv_v));
	memset(&curr, 0, sizeof(fmintv_v));
	memset(&tmp, 0, sizeof(fmintv6_v));
	for (i = w->start; i < w->n; i += w->step) {
		int len = w->len[i];
		uint8_t *seq = w->seq[i], *qual = w->qual[i];
		for (j = 0; j < len; ++j) {
			seq[j] = seq[j] < 128? seq_nt6_table[seq[j]] : 5;
			qual[j] -= 33;
		}
		fm6_ec2_corr1(w->opt, w->e, len, seq, qual, &prev, &curr, &tmp);
		for (j = 0; j < len; ++j) {
			seq[j] = "$ACGTN"[seq[j]];
			qual[j] += 33;
		}
	}
	free(prev.a); free(curr.a); free(tmp.a);
	return 0;
}

static void batch_process(const fmcopt_t *opt, const rld_t *e, int64_t tot, int n, const int *len, uint8_t **seq, uint8_t **qual)
{
	pthread_t *tid;
	worker_t *w;
	int j;
	kstring_t str;

	w = (worker_t*)calloc(opt->n_threads, sizeof(worker_t));
	tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
	for (j = 0; j < opt->n_threads; ++j) {
		worker_t *ww = w + j;
		ww->e = e;
		ww->opt = opt;
		ww->step = opt->n_threads;
		ww->start = j;
		ww->len = len; ww->seq = seq; ww->qual = qual;
	}
	for (j = 0; j < opt->n_threads; ++j) pthread_create(&tid[j], 0, worker, w + j);
	for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
	str.l = str.m = 0; str.s = 0;
	for (j = 0; j < n; ++j) {
		str.l = 0;
		kputc('@', &str); kputl(tot + j, &str); kputc('\n', &str);
		kputs((char*)seq[j], &str); kputsn("\n+\n", 3, &str);
		kputs((char*)qual[j], &str);
		puts(str.s);
	}
	free(str.s);
	free(w); free(tid);
}

int fm6_correct2(const fmcopt_t *opt, const rld_t *e, const char *fn)
{
	kseq_t *kseq;
	gzFile fp;
	int64_t tot_len, tot;
	int n, m, *len;
	uint8_t **seq, **qual;

	fp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	kseq = kseq_init(fp);
	n = m = 0; tot_len = 0; tot = 0;
	seq = qual = 0; len = 0;
	while (kseq_read(kseq) >= 0) {
		if (tot_len + kseq->seq.l > opt->chunk_size) {
			batch_process(opt, e, tot, n, len, seq, qual);
			n = 0;
		}
		if (n == m) {
			m = m? m<<1 : 256;
			len  = realloc(len,  m * sizeof(int));
			seq  = realloc(seq,  m * sizeof(void*));
			qual = realloc(qual, m * sizeof(void*));
		}
		len[n]  = kseq->seq.l;
		seq[n]  = (uint8_t*)strdup(kseq->seq.s);
		qual[n] = (uint8_t*)strdup(kseq->qual.s);
		++n; tot_len += kseq->seq.l; ++tot;
	}
	batch_process(opt, e, tot, n, len, seq, qual);
	free(len); free(seq); free(qual);
	kseq_destroy(kseq);
	gzclose(fp);
	return 0;
}

#include <unistd.h>

int main_fmc(int argc, char *argv[])
{
	int c, use_mmap = 0;
	fmcopt_t opt;
	rld_t *e;
	fmc_opt_init(&opt);
	while ((c = getopt(argc, argv, "Ml:t:")) >= 0) {
		switch (c) {
			case 'l': opt.min_l = atoi(optarg); break;
			case 't': opt.n_threads = atoi(optarg); break;
			case 'M': use_mmap = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi fmc [options] <idx.fmd> <reads.fq>\n\n");
		fprintf(stderr, "Options: -l INT      min match [%d]\n", opt.min_l);
		fprintf(stderr, "         -t INT      number of threads [1]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	fm6_correct2(&opt, e, argv[optind+1]);
	rld_destroy(e);
	return 0;
}
