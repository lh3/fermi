#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <limits.h>
#include <ctype.h>
#include "rld.h"
#include "fermi.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int sais(const unsigned char *T, int *SA, int n, int k);
int sais64(const unsigned char *T, int64_t *SA, int64_t n, int64_t k);
int sa_check(const unsigned char *T, const int *SA, int n);

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
	int i, j, l, c = 0, plain = 0, use_mmap = 0;
	uint64_t *cnt, *rank, sum = 0;
	double t;
	while ((c = getopt(argc, argv, "PM")) >= 0) {
		switch (c) {
			case 'P': plain = 1; break;
			case 'M': use_mmap = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fermi chkbwt [-MP] <idxbase.bwt>\n");
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
		if (plain) for (i = 0; i < l; ++i) putchar("$ACGTN"[c]);
	}
	fprintf(stderr, "[M::%s] Checked the rank function in %.3lf seconds.\n", __func__, cputime() - t);
	for (j = 0; j < e->asize; ++j) {
		if (cnt[j] != e->mcnt[j+1]) {
			fprintf(stderr, "[E::%s] Different counts for symbol %d: %lld != %lld\n", __func__, j,
					(unsigned long long)cnt[j], (unsigned long long)e->mcnt[j+1]);
		}
	}
	if (plain) putchar('\n');
	rld_destroy(e);
	return 0;
}

static void print_i(const rld_t *e, uint64_t i, kstring_t *s)
{
	int j;
	fm_retrieve(e, i, s);
	for (j = s->l - 1; j >= 0; --j)
		putchar("$ACGTN"[(int)s->s[j]]);
	putchar('\n');
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
		fprintf(stderr, "Usage: fermi unpack [-M] [-i index] <idxbase.bwt>\n");
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

int main_exact(int argc, char *argv[])
{
	int c, min_match = 30, use_mmap = 0;
	rld_t *e;
	kseq_t *seq;
	gzFile fp;
	while ((c = getopt(argc, argv, "m:M")) >= 0) {
		switch (c) {
			case 'm': min_match = atoi(optarg); break;
			case 'M': use_mmap = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fermi exact [-M] [-m minMatch] <idxbase.bwt> <src.fa>\n");
		return 1;
	}
	fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	while (kseq_read(seq) >= 0) {
		seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
		printf(">%s\n", seq->name.s);
		printf("%d\n", fm6_search_overlap(e, min_match, seq->seq.l, (uint8_t*)seq->seq.s, 0));
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
		fprintf(stderr, "Usage: fermi merge [-fh] [-o out.bwt] [-t nThreads] <in0.bwt> <in1.bwt> [...]\n");
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
	int sbits = 3, plain = 0, force = 0, no_reverse = 0, asize = 6, inc_N = 0;
	int64_t i, sum_l = 0, l, max, block_size = 250000000;
	uint8_t *s;
	char *idxfn = 0;
	double t;
	rld_t *e = 0;

	{ // parse the command line
		int c;
		while ((c = getopt(argc, argv, "NPRfdb:o:i:s:")) >= 0) {
			switch (c) {
				case 'i':
					e = rld_restore(optarg);
					if (e == 0) {
						fprintf(stderr, "[E::%s] Fail to open the index file `%s'.\n", __func__, optarg);
						return 1;
					}
					break;
				case 'P': plain = 1; break;
				case 'f': force = 1; break;
				case 'R': no_reverse = 1; break;
				case 'b': sbits = atoi(optarg); break;
				case 'o': idxfn = strdup(optarg); break;
				case 's': block_size = atol(optarg); break;
				case 'N': inc_N = 1; break;
			}
		}
		if (argc == optind) {
			fprintf(stderr, "Usage: fermi build [-fRP] [-b sbits] [-o outFile] <in.fa>\n");
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
		gzFile fp;
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
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
			if (!no_reverse) {
				seq_revcomp6(seq->seq.l, (uint8_t*)seq->seq.s);
				memcpy(s + l, seq->seq.s, seq->seq.l + 1);
				l += seq->seq.l + 1;
			}
			sum_l += (seq->seq.l + 1) * 2;
		}
		kseq_destroy(seq);
		gzclose(fp);
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
