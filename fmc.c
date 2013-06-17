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
	int max_pen, max_d, diff_factor;
	double ratio_factor;
} fmcopt_t;

typedef struct {
	fmintv_t k[6];
} fmintv6_t;

typedef struct { size_t n, m; fmintv6_t *a; } fmintv6_v;

void fmc_opt_init(fmcopt_t *opt)
{
	memset(opt, 0, sizeof(fmcopt_t));
	opt->chunk_size = 1<<26;
	opt->n_threads = 1;
	opt->min_l = 23;
	opt->min_occ = 6;
	opt->for_qthres = 30;
	opt->rev_qthres = 20;
	opt->multi_thres = 10;

	opt->max_pen = 60;
	opt->max_d = 7;
	opt->diff_factor = 13;
	opt->ratio_factor = 10.;
}

static inline void ec_cns_gen(const fmcopt_t *opt, const fmintv_t k[6], uint8_t q[4])
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

static inline void ec_cns_upd(uint8_t q[4], const uint8_t q1[4], int is_comp)
{
	int c;
	if (is_comp) {
		for (c = 0; c < 4; ++c)
			q[c] = q[c] < q1[3-c]? q[c] : q1[3-c];
	} else {
		for (c = 0; c < 4; ++c)
			q[c] = q[c] < q1[c]? q[c] : q1[c];
	}
}

static inline int ec_cns_call(const uint8_t q[4], int base, int qual)
{
	int min, min_c, min2, c, iq[4];
	for (c = 0; c < 4; ++c) iq[c] = q[c];
	if (base >= 0 && base <= 3) iq[base] -= qual;
	for (c = 0, min_c = -1, min = min2 = 255; c < 4; ++c)
		if (min > iq[c]) min2 = min, min = iq[c], min_c = c;
		else if (min2 > iq[c]) min2 = iq[c];
	return min_c << 8 | (min2 - min);
}

typedef struct {
	uint8_t pen[4];
	uint8_t ob, oq, eb, eq;
} ecseq_t;

static ecseq_t *ec_seq_gen(int l_seq, const uint8_t *seq, const uint8_t *qual)
{
	int i;
	ecseq_t *s;
	s = calloc(l_seq, sizeof(ecseq_t));
	for (i = 0; i < l_seq; ++i) {
		uint32_t *p = (uint32_t*)s[i].pen;
		*p = 0xffffffffU;
		s[i].ob = seq[i]; s[i].eq = s[i].oq = qual[i];
	}
	return s;
}

static void ec_seq_rev(int l_seq, ecseq_t *seq)
{
	int i;
	for (i = 0; i < l_seq>>1; ++i) {
		ecseq_t tmp, *s = &seq[i], *t = &seq[l_seq - 1 - i];
		tmp = *s; *s = *t; *t = tmp;
	}
	for (i = 0; i < l_seq; ++i) {
		ecseq_t *si = &seq[i];
		uint32_t *p = (uint32_t*)si->pen;
		si->ob = fm6_comp(si->ob); si->eb = fm6_comp(si->eb);
		*p = (*p & 0x0000FFFFU) << 16 | *p >> 16;
		*p = (*p & 0x00FF00FFU) << 8  | (*p & 0xFF00FF00U) >> 8;
	}
}

int fm6_ec2_core(const fmcopt_t *opt, const rld_t *e, int l_seq, ecseq_t *seq, int x, fmintv_v *prev, fmintv_v *curr, fmintv6_v *tmp)
{ // this function is similar to fm6_smem1_core()
	int i, j, ret;
	fmintv_t ik;
	fmintv_v *swap;
	fmintv6_t ok;

	//fprintf(stderr, "x=%d\n", x);
	fm6_set_intv(e, seq[x].ob, ik);
	ik.info = x + 1;
	for (i = x + 1, curr->n = 0; i < l_seq; ++i) { // forward search
		ecseq_t *si = &seq[i];
		int c;
		fm6_extend(e, &ik, ok.k, 0); // forward extension
		if (i - x >= opt->min_l) { // then check the multi-alignment
			int call;
			uint8_t q[4];
			ec_cns_gen(opt, ok.k, q);
			ec_cns_upd(si->pen, q, 1);
			call = ec_cns_call(si->pen, fm6_comp(si->ob) - 1, si->oq);
			//fprintf(stderr, "[F,%d,%c]\tpen[4]={%d,%d,%d,%d}\teb=%c\teq=%d\n", i, "$ACGTN"[si->ob], q[0], q[1], q[2], q[3], "ACGT"[call>>8], call&0xff);
			si->eb = (call >> 8) + 1;
			si->eq = call & 0xff;
		}
		c = si->eb? si->eb : si->ob;
		c = fm6_comp(c);
		if (ok.k[c].x[2] != ik.x[2]) { // change of interval size
			kv_push(fmintv_t, *curr, ik);
			if (ok.k[c].x[2] < opt->min_occ) break;
		}
		ik = ok.k[c]; ik.info = i + 1;
	}

	if (i == l_seq) kv_push(fmintv_t, *curr, ik);
	fm_reverse_fmivec(curr);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= 0; --i) { // backward search
		ecseq_t *si = &seq[i];
		int c, call, is_ec;
		kv_resize(fmintv6_t, *tmp, prev->n);
		tmp->n = prev->n;
		for (j = 0; j < prev->n; ++j) // collect all the intervals at the current position
			fm6_extend(e, &prev->a[j], tmp->a[j].k, 1);
		for (is_ec = 0, j = 0; j < prev->n; ++j) { // check if we need to make a correction
			if (prev->a[j].info - i >= opt->min_l) {
				uint8_t q[4];
				ec_cns_gen(opt, tmp->a[j].k, q);
				ec_cns_upd(si->pen, q, 0);
				is_ec = 1;
			}
		}
		if (is_ec) {
			call = ec_cns_call(si->pen, si->ob - 1, si->oq);
			si->eb = (call >> 8) + 1;
			si->eq = call & 0xff;
		}
		c = si->eb? si->eb : si->ob;
		for (j = 0, curr->n = 0; j < prev->n; ++j) { // update $curr
			fmintv_t *ok = tmp->a[j].k;
			if (ok[c].x[2] >= opt->min_occ && (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2])) {
				ok[c].info = prev->a[j].info;
				kv_push(fmintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	//fprintf(stderr, "ret=%d\n", ret);

	return ret;
}

static int fm6_ec2_corr1(const fmcopt_t *opt, const rld_t *e, int l_seq, uint8_t *seq, uint8_t *qual, fmintv_v *prev, fmintv_v *curr, fmintv6_v *tmp)
{
	int x;
	ecseq_t *s;
	s = ec_seq_gen(l_seq, seq, qual);
	// error correction
	x = 0;
	while ((x = fm6_ec2_core(opt, e, l_seq, s, x, prev, curr, tmp)) < l_seq);
	ec_seq_rev(l_seq, s);
	x = 0;
	while ((x = fm6_ec2_core(opt, e, l_seq, s, x, prev, curr, tmp)) < l_seq);
	ec_seq_rev(l_seq, s);
	// modify $seq and $qual
	for (x = 0; x < l_seq; ++x) {
		ecseq_t *sx = &s[x];
		seq[x] = sx->eb == 0 || sx->ob == sx->eb? "$ACGTN"[sx->ob] : "$acgtn"[sx->eb];
		qual[x] = sx->eq + 33;
	}
	free(s);
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
		ww->n = n; ww->len = len; ww->seq = seq; ww->qual = qual;
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
