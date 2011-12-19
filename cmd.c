#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <limits.h>
#include <ctype.h>
#include "rld.h"
#include "fermi.h"
#include "kstring.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int sais(const unsigned char *T, int *SA, int n, int k);
int sais64(const unsigned char *T, int64_t *SA, int64_t n, int64_t k);

void seq_char2nt6(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

double cputime();

int main_cnt2qual(int argc, char *argv[])
{
	int q = 17, i;
	gzFile fp;
	kseq_t *seq;
	if (argc < 2) {
		fprintf(stderr, "Usage: fermi cnt2qual <in.fq> [%d]\n", q);
		return 1;
	}
	if (argc >= 3) q = atoi(argv[2]);
	fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		if (seq->qual.l) {
			for (i = 0; i < seq->qual.l; ++i) {
				int x = q * (seq->qual.s[i] - 33) + 33;
				seq->qual.s[i] = x > 126? 126 : x;
			}
		}
		putchar('@'); fputs(seq->name.s, stdout);
		if (seq->comment.l) {
			putchar('\t'); puts(seq->comment.s);
		} else putchar('\n');
		puts(seq->seq.s);
		if (seq->qual.l) {
			putchar('+'); putchar('\n');
			puts(seq->qual.s);
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int main_chkbwt(int argc, char *argv[])
{
	rld_t *e;
	rlditr_t itr;
	int i, j, l, c = 0, plain = 0, use_mmap = 0, check_rank = 0;
	uint64_t *cnt, *rank, sum = 0;
	double t;
	while ((c = getopt(argc, argv, "pMr")) >= 0) {
		switch (c) {
			case 'p': plain = 1; break;
			case 'r': check_rank = 1; break;
			case 'M': use_mmap = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi chkbwt [options] <idxbase.bwt>\n\n");
		fprintf(stderr, "Options: -M        load the FM-index as a memory mapped file\n");
		fprintf(stderr, "         -r        check rank\n");
		fprintf(stderr, "         -p        print the BWT to the stdout\n\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	if (e == 0) {
		fprintf(stderr, "[E::%s] Fail to read the index file.\n", __func__);
		return 1;
	}
	if (fm_verbose >= 3) {
		fprintf(stderr, "[M::%s] marginal counts: (", __func__);
		for (i = 0; i < e->asize1; ++i)
			fprintf(stderr, "%lld%s", (long long)e->mcnt[i], i == e->asize? ")" : ", ");
		fputc('\n', stderr);
	}
	cnt = alloca(e->asize * 8);
	rank = alloca(e->asize * 8);
	rld_itr_init(e, &itr, 0);
	for (i = 0; i < e->asize; ++i) cnt[i] = 0;
	t = cputime();
	while ((l = rld_dec(e, &itr, &c, 0)) >= 0) {
		if (c >= e->asize) {
			fprintf(stderr, "[E::%s] Symbol `%d' is not in the alphabet.\n", __func__, c);
			exit(1); // memory leak
		}
		if (check_rank) {
			for (i = 0; i < l; ++i) {
				++cnt[c];
				rld_rank1a(e, sum, rank);
				for (j = 0; j < e->asize; ++j) {
					if (cnt[j] != rank[j]) {
						fprintf(stderr, "[E::%s] rank(%d,%lld)=%lld != %lld\n", __func__,
								j, (unsigned long long)sum, (unsigned long long)rank[j], (unsigned long long)cnt[j]);
						exit(1); // memory leak
					}
				}
				++sum;
				if (sum%10000000 == 0)
					fprintf(stderr, "[M::%s] Checked %lld symbols.\n", __func__, (unsigned long long)sum);
			}
		}
		if (plain) for (i = 0; i < l; ++i) putchar("$ACGTN"[c]);
	}
	if (check_rank) {
		fprintf(stderr, "[M::%s] Checked the rank function in %.3lf seconds.\n", __func__, cputime() - t);
		for (j = 0; j < e->asize; ++j) {
			if (cnt[j] != e->mcnt[j+1]) {
				fprintf(stderr, "[E::%s] Different counts for symbol %d: %lld != %lld\n", __func__, j,
						(unsigned long long)cnt[j], (unsigned long long)e->mcnt[j+1]);
			}
		}
	}
	if (plain) putchar('\n');
	rld_destroy(e);
	return 0;
}

static void print_i(const rld_t *e, uint64_t i, kstring_t *s)
{
	int j;
	uint64_t k;
	k = fm_retrieve(e, i, s);
	for (j = s->l - 1; j >= 0; --j)
		putchar("$ACGTN"[(int)s->s[j]]);
	printf("\t%ld\n", (long)k);
}

int main_unpack(int argc, char *argv[])
{
	rld_t *e;
	int c, n, m, use_mmap = 0;
	uint64_t i, *list;
	kstring_t s;
	s.m = s.l = 0; s.s = 0;
	n = m = 0; list = 0;
	while ((c = getopt(argc, argv, "Mi:")) >= 0) {
		switch (c) {
			case 'i':
				if (n == m) {
					m = m? m<<1 : 16;
					list = realloc(list, 8 * m);
				}
				list[n++] = atol(optarg); break;
				break;
			case 'M': use_mmap = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi unpack [-M] [-i index] <seqs.bwt>\n\n");
		fprintf(stderr, "Options: -i INT    index of the read to output, starting from 0 [null]\n");
		fprintf(stderr, "         -M        load the FM-index as a memory mapped file\n\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	if (n) {
		for (i = 0; (int)i < n; ++i)
			if (list[i] >= 0 && list[i] < e->mcnt[1])
				print_i(e, list[i], &s);
	} else {
		for (i = 0; i < e->mcnt[1]; ++i)
			print_i(e, i, &s);
	}
	rld_destroy(e);
	free(s.s);
	return 0;
}

int main_unitig(int argc, char *argv[])
{
	int c, use_mmap = 0, n_threads = 1, min_match = 30;
	rld_t *e;
	uint64_t *sorted = 0;
	char *fn_sorted = 0;
	while ((c = getopt(argc, argv, "Ml:t:r:")) >= 0) {
		switch (c) {
			case 'l': min_match = atoi(optarg); break;
			case 'M': use_mmap = 1; break;
			case 't': n_threads = atoi(optarg); break;
			case 'r': fn_sorted = strdup(optarg); break;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi unitig [options] <reads.bwt>\n\n");
		fprintf(stderr, "Options: -l INT      min match [%d]\n", min_match);
		fprintf(stderr, "         -t INT      number of threads [1]\n");
		fprintf(stderr, "         -r FILE     rank file [null]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	if (fn_sorted) {
		FILE *fp;
		fp = fopen(fn_sorted, "rb");
		sorted = malloc(e->mcnt[1] * 8);
		fread(sorted, e->mcnt[1], 8, fp);
		fclose(fp);
		free(fn_sorted);
	}
	fm6_unitig(e, min_match, n_threads, sorted);
	free(sorted);
	rld_destroy(e);
	return 0;
}

int main_correct(int argc, char *argv[])
{
	int c, use_mmap = 0, n_threads = 1;
	rld_t *e;
	fmecopt_t opt;
	opt.w = 23; opt.min_occ = 3; opt.keep_bad = 0; opt.is_paired = 0; opt.max_corr = 0.3;
	while ((c = getopt(argc, argv, "MKt:k:v:O:pC:")) >= 0) {
		switch (c) {
			case 'M': use_mmap = 1; break;
			case 'K': opt.keep_bad = 1; break;
			case 't': n_threads = atoi(optarg); break;
			case 'k': opt.w = atoi(optarg); break;
			case 'v': fm_verbose = atoi(optarg); break;
			case 'O': opt.min_occ = atoi(optarg); break;
			case 'p': opt.is_paired = 1; break;
			case 'C': opt.max_corr = atof(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi correct [options] <reads.fmd> <reads.fq>\n\n");
		fprintf(stderr, "Options: -k INT      k-mer length [%d]\n", opt.w);
		fprintf(stderr, "         -O INT      minimum (k+1)-mer occurrences [%d]\n", opt.min_occ);
		fprintf(stderr, "         -t INT      number of threads [%d]\n", n_threads);
		fprintf(stderr, "         -C FLOAT    max fraction of corrected bases [%.2f]\n", opt.max_corr);
		fprintf(stderr, "         -K          keep bad/unfixable reads\n");
		fprintf(stderr, "\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	fm6_ec_correct(e, &opt, argv[optind+1], n_threads);
	rld_destroy(e);
	return 0;
}

int main_exact(int argc, char *argv[])
{
	int c, i, use_mmap = 0, self_match = 0;
	rld_t *e;
	kseq_t *seq;
	gzFile fp;
	kstring_t str;
	fmintv_v a;

	while ((c = getopt(argc, argv, "Ms")) >= 0) {
		switch (c) {
			case 'M': use_mmap = 1; break;
			case 's': self_match = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fermi exact [-Ms] <idxbase.bwt> <src.fa>\n");
		return 1;
	}
	fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);

	a.m = a.n = 0; a.a = 0;
	str.m = str.l = 0; str.s = 0;
	while (kseq_read(seq) >= 0) {
		seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
		fm6_smem(e, seq->seq.l, (uint8_t*)seq->seq.s, &a, self_match);
		str.l = 0; kputs("SQ\t", &str); kputs(seq->name.s, &str); kputc('\t', &str); kputw(seq->seq.l, &str); kputc('\t', &str); kputw(a.n, &str);
		puts(str.s);
		for (i = 0; i < a.n; ++i) {
			fputs("EM\t", stdout);
			fm6_write_smem(e, &a.a[i], &str);
			puts(str.s);
		}
		puts("//");
	}
	rld_destroy(e);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int main_merge(int argc, char *argv[])
{
	int i, c, force = 0, n_threads = 1;
	rld_t *e0, *e1;
	char *idxfn = 0;
	while ((c = getopt(argc, argv, "fo:t:")) >= 0) {
		switch (c) {
			case 'f': force = 1; break;
			case 't': n_threads = atoi(optarg); break;
			case 'o': idxfn = strdup(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi merge [-f] [-o out.bwt] [-t nThreads] <in0.bwt> <in1.bwt> [...]\n\n");
		fprintf(stderr, "Options: -f        force to overwrite the output file (effective with -o)\n");
		fprintf(stderr, "         -o FILE   output file name [null]\n");
		fprintf(stderr, "         -t INT    number of threads to use\n\n");
		return 1;
	}
	if (force == 0 && idxfn) {
		FILE *fp = fopen(idxfn, "r");
		if (fp) {
			fclose(fp);
			fprintf(stderr, "[E::%s] File `%s' exists. Please use `-f' to overwrite.\n", __func__, idxfn);
			return 1;
		}
	}
	if (idxfn == 0) idxfn = strdup("-");
	e0 = rld_restore(argv[optind]);
	fprintf(stderr, "[M::%s] Loaded file `%s'.\n", __func__, argv[optind]);
	for (i = optind + 1; i < argc; ++i) {
		e1 = rld_restore(argv[i]);
		fprintf(stderr, "[M::%s] Loaded file `%s'.\n", __func__, argv[i]);
		e0 = fm_merge(e0, e1, n_threads); // e0 and e1 will be deallocated during merge
		fprintf(stderr, "[M::%s] Merged file `%s' to the existing index.\n", __func__, argv[i]);
	}
	rld_dump(e0, idxfn);
	rld_destroy(e0);
	free(idxfn);
	return 0;
}

int main_build(int argc, char *argv[]) // this routinue to replace main_index() in future
{
	int sbits = 3, force = 0, asize = 6;
	int64_t sum_l = 0, l, max, block_size = 250000000;
	uint8_t *s;
	char *idxfn = 0;
	double t;
	rld_t *e = 0;

	{ // parse the command line
		int c;
		while ((c = getopt(argc, argv, "fb:o:i:s:")) >= 0) {
			switch (c) {
				case 'i':
					e = rld_restore(optarg);
					if (e == 0) {
						fprintf(stderr, "[E::%s] Fail to open the index file `%s'.\n", __func__, optarg);
						return 1;
					}
					break;
				case 'f': force = 1; break;
				case 'b': sbits = atoi(optarg); break;
				case 'o': idxfn = strdup(optarg); break;
				case 's': block_size = atol(optarg); break;
			}
		}
		if (argc == optind) {
			fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   fermi build [options] <in.fa>\n\n");
			fprintf(stderr, "Options: -b INT    use a small marker per 2^(INT+3) bytes [%d]\n", sbits);
			fprintf(stderr, "         -f        force to overwrite the output file (effective with -o)\n");
			fprintf(stderr, "         -i FILE   append the FM-index to the existing FILE [null]\n");
			fprintf(stderr, "         -o FILE   output file name [null]\n");
			fprintf(stderr, "         -s INT    number of symbols to process at a time [%ld]\n", (long)block_size);
			fprintf(stderr, "\n");
			return 1;
		}
		if (idxfn) {
			FILE *fp;
			if (!force && (fp = fopen(idxfn, "rb")) != 0) {
				fclose(fp);
				fprintf(stderr, "[E::%s] File `%s' exists. Please use `-f' to overwrite.\n", __func__, idxfn);
				return 1;
			}
		} else idxfn = strdup("-");
	}
	
	{ // read sequences
		kseq_t *seq;
		kstring_t str;
		gzFile fp;

		str.l = str.m = 0; str.s = 0;
		fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
		if (fp == 0) {
			fprintf(stderr, "[E::%s] Fail to open the input file.\n", __func__);
			return 1;
		}
		seq = kseq_init(fp);
		l = 0; max = 16;
		s = malloc(max);
		t = cputime();
		while (kseq_read(seq) >= 0) {
			if (l && l + (seq->seq.l + 1) * 2 > block_size) {
				e = fm_build(e, asize, sbits, l, s);
				fprintf(stderr, "[M::%s] Constructed BWT for %lld million symbols in %.3f seconds.\n", __func__, (long long)sum_l/1000000, cputime() - t);
				l = 0;
			}
			if (l + (seq->seq.l + 1) * 2 > max) { // we do not set max as block_size because this is more flexible
				max = l + (seq->seq.l + 1) * 2 + 1;
				kroundup32(max);
				s = realloc(s, max);
			}
			seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
			seq_revcomp6(seq->seq.l, (uint8_t*)seq->seq.s);
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
			sum_l += (seq->seq.l + 1) * 2;
		}
		kseq_destroy(seq);
		gzclose(fp);
		free(str.s);
		if (l) {
			e = fm_build(e, asize, sbits, l, s);
			fprintf(stderr, "[M::%s] Constructed BWT for %lld million symbols in %.3f seconds.\n", __func__, (long long)sum_l/1000000, cputime() - t);
		}
	}

	rld_dump(e, idxfn);
	rld_destroy(e);
	free(s); free(idxfn);
	return 0;
}

int main_seqsort(int argc, char *argv[])
{
	extern uint64_t *fm6_seqsort(const rld_t *e, int n_threads);
	int c, n_threads = 1;
	rld_t *e;
	uint64_t *sorted;
	while ((c = getopt(argc, argv, "t:")) >= 0) {
		switch (c) {
			case 't': n_threads = atoi(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi seqsort [-t nThreads=1] <reads.fmd>\n");
		return 1;
	}
	e = rld_restore(argv[optind]);
	sorted = fm6_seqsort(e, n_threads);
	fwrite(sorted, 8, e->mcnt[1], stdout);
	free(sorted);
	rld_destroy(e);
	return 0;
}

int main_peread(int argc, char *argv[])
{
	extern int msg_peread(const msg_t *g, double avg, double std);
	int c;
	double avg, std;
	msg_t *g;

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 3 > argc) {
		fprintf(stderr, "Usage: fermi peread <in.fmd> <avg> <std>\n");
		return 1;
	}
	avg = atof(argv[optind+1]);
	std = atof(argv[optind+2]);
	g = msg_read(argv[optind], 0, 1<<30, 0.);
	msg_peread(g, avg, std);
	msg_destroy(g);
	return 0;
}

int main_clean(int argc, char *argv[])
{
	msg_t *g;
	int c, do_clean = 0, max_arc = 512;
	float read_diff_ratio = 0.7;
	fmclnopt_t opt;
	opt.aggressive_pop = 0;
	opt.min_tip_len = 200;
	opt.min_weak_cov= 1.01;
	opt.min_bub_cov = 10.; opt.min_bub_ratio= 0.3;
	opt.min_ovlp    = 60;  opt.min_ovlp_ratio=0.8;
	opt.n_iter = 3;
	while ((c = getopt(argc, argv, "CAl:c:T:r:w:o:R:n:N:d:")) >= 0) {
		switch (c) {
			case 'l': opt.min_tip_len =  atoi(optarg); break;
			case 'c': opt.min_weak_cov=  atof(optarg); break;
			case 'w': opt.min_bub_cov =  atof(optarg); break;
			case 'r': opt.min_bub_ratio= atof(optarg); break;
			case 'o': opt.min_ovlp    =  atoi(optarg); break;
			case 'R': opt.min_ovlp_ratio=atof(optarg); break;
			case 'n': opt.n_iter = atoi(optarg); break;
			case 'N': max_arc = atoi(optarg); break;
			case 'd': read_diff_ratio =  atof(optarg); break;
			case 'C': do_clean = 1; break;
			case 'A': opt.aggressive_pop = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi clean [options] <in.msg>\n\n");
		fprintf(stderr, "Options: -N INT      read maximum INT neighbors per node [%d]\n", max_arc);
		fprintf(stderr, "         -d FLOAT    drop a neighbor if relative overlap ratio below FLOAT [%.2f]\n\n", read_diff_ratio); 
		fprintf(stderr, "         -C          clean the graph\n");
		fprintf(stderr, "         -l INT      minimum tip length [%d]\n", opt.min_tip_len);
		fprintf(stderr, "         -o INT      minimum overlap [%d]\n", opt.min_ovlp);
		fprintf(stderr, "         -R FLOAT    minimum relative overlap ratio [%.2f]\n", opt.min_ovlp_ratio);
		fprintf(stderr, "         -c FLOAT    minimum node coverage [%.1f]\n", opt.min_weak_cov);
		fprintf(stderr, "         -n INT      number of iterations [%d]\n", opt.n_iter);
		fprintf(stderr, "         -A          aggressive bubble popping\n");
		fprintf(stderr, "         -w FLOAT    minimum simple bubble coverage [%.1f]\n", opt.min_bub_cov);
		fprintf(stderr, "         -r FLOAT    minimum simple bubble ratio [%.2f]\n", opt.min_bub_ratio);
		fprintf(stderr, "\n");
		return 1;
	}
	g = msg_read(argv[optind], 1, max_arc, read_diff_ratio);
	msg_join_unambi(g);
	if (do_clean) msg_clean(g, &opt);
	msg_print(&g->nodes);
	msg_destroy(g);
	return 0;
}
