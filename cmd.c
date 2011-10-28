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

int main_strlen(int argc, char *argv[])
{
	int64_t l = 0;
	gzFile fp;
	kseq_t *seq;
	if (argc == 1) {
		fprintf(stderr, "Usage: fermi strlen <in.fq>\n");
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) l += 2 * (seq->seq.l + 1);
	printf("%lld\n", (long long)l);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int main_chkbwt(int argc, char *argv[])
{
	rld_t *e;
	rlditr_t itr;
	int i, j, l, c = 0, plain = 1, use_mmap = 0, check_rank = 0;
	uint64_t *cnt, *rank, sum = 0;
	double t;
	while ((c = getopt(argc, argv, "PMr")) >= 0) {
		switch (c) {
			case 'P': plain = 0; break;
			case 'r': check_rank = 1; break;
			case 'M': use_mmap = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi chkbwt [-MPrR] <idxbase.bwt>\n\n");
		fprintf(stderr, "Options: -M        load the FM-index as a memory mapped file\n");
		fprintf(stderr, "         -r        check rank\n");
		fprintf(stderr, "         -P        do not print the BWT to the stdout\n\n");
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
	while ((c = getopt(argc, argv, "SMDl:t:r:m:L:")) >= 0) {
		switch (c) {
			case 'l': min_match = atoi(optarg); break;
			case 'M': use_mmap = 1; break;
			case 't': n_threads = atoi(optarg); break;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi unitig [options] <reads.bwt>\n\n");
		fprintf(stderr, "Options: -l INT      min match [%d]\n", min_match);
		fprintf(stderr, "         -t INT      number of threads [1]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	fm6_unitig(e, min_match, n_threads);
	rld_destroy(e);
	return 0;
}

int main_correct(int argc, char *argv[])
{
	int c, use_mmap = 0, n_threads = 1, _w, _T;
	rld_t *e;
	fmecopt_t opt;
	opt.cov = 30.0; opt.t = 3; opt.T = opt.w = 0; opt.err = 0.01; opt.max_pre_mm = 8;
	while ((c = getopt(argc, argv, "A:Mt:k:T:c:m:e:v:")) >= 0) {
		switch (c) {
			case 'M': use_mmap = 1; break;
			case 'm': opt.t = atoi(optarg); break;
			case 'c': opt.cov = atof(optarg); break;
			case 'e': opt.err = atof(optarg); break;
			case 't': n_threads = atoi(optarg); break;
			case 'k': opt.w = atoi(optarg); break;
			case 'T': opt.T = atoi(optarg); break;
			case 'v': fm_verbose = atoi(optarg); break;
			case 'A': opt.max_pre_mm = atoi(optarg); break;
		}
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi correct [options] <reads.bwt>\n\n");
		fprintf(stderr, "Options: -c FLOAT    expected coverage [%.1f]\n", opt.cov);
		fprintf(stderr, "         -e FLOAT    expected per-base error rate [%.2f]\n", opt.err);
		fprintf(stderr, "         -m INT      do not correct an error appearing more than INT times [%d]\n", opt.t);
		fprintf(stderr, "         -k INT      k-mer length [inferred from -c/-e]\n");
		fprintf(stderr, "         -T INT      threshold for a correct base [inferred from -c/-e]\n");
		fprintf(stderr, "         -t INT      number of threads [%d]\n", n_threads);
		fprintf(stderr, "         -A INT      max #clustered errors in prefix [%d]\n\n", opt.max_pre_mm);
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	fm_ec_genpar(e->mcnt[1]/2, (int)((double)e->mcnt[0] / e->mcnt[1] - 1 + 0.5), opt.cov, opt.err, &_w, &_T);
	if (opt.w <= 0) opt.w = _w;
	if (opt.T <= 0) opt.T = _T;
	fm6_ec_correct(e, &opt, n_threads);
	rld_destroy(e);
	return 0;
}

int main_exact(int argc, char *argv[])
{
	int c, i, use_mmap = 0;
	rld_t *e;
	kseq_t *seq;
	gzFile fp;
	kstring_t str;
	fmintv_v a;

	while ((c = getopt(argc, argv, "M")) >= 0) {
		switch (c) {
			case 'M': use_mmap = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fermi exact [-M] <idxbase.bwt> <src.fa>\n");
		return 1;
	}
	fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);

	a.m = a.n = 0; a.a = 0;
	str.m = str.l = 0; str.s = 0;
	while (kseq_read(seq) >= 0) {
		seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
		fm6_smem(e, seq->seq.l, (uint8_t*)seq->seq.s, &a);
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
	int sbits = 3, plain = 0, force = 0, asize = 6, inc_N = 0, min_q = 0, trim_end_N = 1;
	int64_t i, sum_l = 0, l, max, block_size = 250000000;
	uint8_t *s;
	char *idxfn = 0;
	double t;
	rld_t *e = 0;

	{ // parse the command line
		int c;
		while ((c = getopt(argc, argv, "NPfb:o:i:s:q:T")) >= 0) {
			switch (c) {
				case 'i':
					e = rld_restore(optarg);
					if (e == 0) {
						fprintf(stderr, "[E::%s] Fail to open the index file `%s'.\n", __func__, optarg);
						return 1;
					}
					break;
				case 'q': min_q = atoi(optarg); break;
				case 'P': plain = 1; break;
				case 'f': force = 1; break;
				case 'b': sbits = atoi(optarg); break;
				case 'o': idxfn = strdup(optarg); break;
				case 's': block_size = atol(optarg); break;
				case 'N': inc_N = 1; break;
				case 'T': trim_end_N = 0; break;
			}
		}
		if (argc == optind) {
			fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   fermi build [-fPN] [-i inBWT] [-b sbits] [-o outFile] [-s blkSize] <in.fa>\n\n");
			fprintf(stderr, "Options: -b INT    use a small marker per 2^(INT+3) bytes [%d]\n", sbits);
			fprintf(stderr, "         -f        force to overwrite the output file (effective with -o)\n");
			fprintf(stderr, "         -i FILE   append the FM-index to the existing FILE [null]\n");
			fprintf(stderr, "         -N        do NOT discard sequences containing 'N' (after trimming)\n");
			fprintf(stderr, "         -o FILE   output file name [null]\n");
			fprintf(stderr, "         -P        output BWT to stdout; do not dump the FM-index\n");
			fprintf(stderr, "         -q INT    convert base with base quality smaller than INT to 'N' [0]\n");
			fprintf(stderr, "         -s INT    number of symbols to process at a time [%ld]\n", (long)block_size);
			fprintf(stderr, "         -T        do NOT trim 'N' at 5'- or 3'-end\n");
			fprintf(stderr, "\n");
			return 1;
		}
		if (!plain) {
			if (idxfn) {
				FILE *fp;
				if (!force && (fp = fopen(idxfn, "rb")) != 0) {
					fclose(fp);
					fprintf(stderr, "[E::%s] File `%s' exists. Please use `-f' to overwrite.\n", __func__, idxfn);
					return 1;
				}
			} else idxfn = strdup("-");
		}
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
			int j, k;
			if (min_q > 0 && seq->qual.l) // convert low-quality bases to 'N'
				for (j = 0; j < seq->seq.l; ++j)
					if (isalpha(seq->seq.s[j]) && seq->qual.s[j] - 33 < min_q)
						seq->seq.s[j] = 'N';
			if (trim_end_N) {
				int end, start, ostart, oend = -1, min_l = seq->seq.l;
				ks_resize(&str, seq->seq.m);
				str.l = 0;
				while (1) {
					for (j = oend + 1; j < seq->seq.l && !isalpha(seq->seq.s[j]); ++j); // skip leading rubbish
					if ((ostart = j) >= seq->seq.l) break; // start of the sequence
					for (j = ostart; toupper(seq->seq.s[j]) == 'N'; ++j); // skip leading N
					start = j;
					for (j = start; isalpha(seq->seq.s[j]); ++j); // find the end
					oend = j;
					for (j = oend - 1; j >= ostart && toupper(seq->seq.s[j]) == 'N'; --j); // skip trailing N
					end = j + 1;
					for (j = start; j < end; ++j) str.s[str.l++] = seq->seq.s[j]; // copy sequence; can be done in place actually
					if (min_l > end - start) min_l = end - start;
					str.s[str.l++] = '.';
					str.s[str.l] = 0;
				}
				if (min_l < 2) continue; // 1bp read? skip
				str.s[str.l-1] = 0;
				memcpy(seq->seq.s, str.s, str.l);
				seq->seq.l = str.l - 1;
			}
			if (!inc_N) { // skip reads containing 'N'
				for (j = 0; j < seq->seq.l; ++j)
					if (toupper(seq->seq.s[j]) == 'N') break;
				if (j != seq->seq.l) continue;
			}
			for (j = k = 0; j < seq->seq.l; ++j) // skip leading non-alphabet characters
				if (isalpha(seq->seq.s[j])) break;
			for (; j < seq->seq.l; ++j) {
				int c = seq->seq.s[j];
				if (!isalpha(c)) c = seq->seq.s[j] = 0;
				if (c || seq->seq.s[k-1]) seq->seq.s[k++] = c;
			}
			seq->seq.l = k; // NB: quality is untouched
			seq->seq.s[k] = 0;
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
			for (j = k = 0; j <= seq->seq.l; ++j) {
				if (seq->seq.s[j] == 0) {
					memcpy(s + l, seq->seq.s + k, j - k + 1);
					l += j - k + 1;
					seq_revcomp6(j - k, (uint8_t*)seq->seq.s + k);
					memcpy(s + l, seq->seq.s + k, j - k + 1);
					l += j - k + 1;
					k = j + 1;
				}
			}
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

	if (plain) {
		for (i = 0; i < l; ++i) putchar("$ACGTN"[s[i]]);
		putchar('\n');
	} else {
		rld_dump(e, idxfn);
		rld_destroy(e);
	}
	free(s); free(idxfn);
	return 0;
}

int main_clean(int argc, char *argv[])
{
	fmnode_v *nodes;
	int c;
	fmclnopt_t opt;
	opt.min_tip_len = 150; opt.min_tip_cov  = 10.0;
	opt.min_bub_cov = 7.0; opt.min_bub_ratio= 0.3;
	opt.min_ovlp    = 0;   opt.min_ovlp_ratio=0.7;
	opt.check = 0;
	while ((c = getopt(argc, argv, "Cl:c:T:r:w:o:R:")) >= 0) {
		switch (c) {
			case 'l': opt.min_tip_len =  atoi(optarg); break;
			case 'c': opt.min_tip_cov =  atof(optarg); break;
			case 'w': opt.min_bub_cov =  atof(optarg); break;
			case 'r': opt.min_bub_ratio= atof(optarg); break;
			case 'o': opt.min_ovlp    =  atoi(optarg); break;
			case 'R': opt.min_ovlp_ratio=atof(optarg); break;
			case 'C': opt.check = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi clean [options] <in.msg>\n\n");
		fprintf(stderr, "Options: -c FLOAT    minimum tip coverage (0 to disable tip trimming) [%.1f]\n", opt.min_tip_cov);
		fprintf(stderr, "         -l INT      minimum tip length [%d]\n", opt.min_tip_len);
		fprintf(stderr, "         -w FLOAT    minimum bubble coverage (0 to disable debubbling) [%.1f]\n", opt.min_bub_cov);
		fprintf(stderr, "         -r FLOAT    minimum bubble ratio [%.2f]\n", opt.min_bub_ratio);
		fprintf(stderr, "         -o INT      minimum overlap [%d]\n", opt.min_ovlp);
		fprintf(stderr, "         -R FLOAT    minimum relative overlap ratio [%.2f]\n", opt.min_ovlp_ratio);
		fprintf(stderr, "\n");
		return 1;
	}
	nodes = msg_read(argv[optind]);
	msg_clean(nodes, &opt);
	msg_print(nodes);
	return 0;
}
