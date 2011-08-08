#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "rld.h"
#include "fermi.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int sais(const unsigned char *T, int *SA, int n, int k);
int sa_check(const unsigned char *T, const int *SA, int n);

void seq_char2nt6(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

double cputime();

int main_index(int argc, char *argv[])
{
	int i, bbits = 3, plain = 0, check = 0, force = 0, no_reverse = 0;
	uint32_t l, max;
	uint8_t *s;
	char *idxfn = 0;

	{ // parse the command line
		int c;
		while ((c = getopt(argc, argv, "CPRfb:")) >= 0) {
			switch (c) {
				case 'C': check = 1; break;
				case 'P': plain = 1; break;
				case 'f': force = 1; break;
				case 'R': no_reverse = 1; break;
				case 'b': bbits = atoi(optarg); break;
			}
		}
		if (argc == optind) {
			fprintf(stderr, "Usage: fermi index [-CfP] [-b sbits] <in.fa>\n");
			return 1;
		}
		if (!plain) {
			FILE *fp;
			idxfn = malloc(strlen(argv[optind]) + 6);
			strcpy(idxfn, argv[optind]);
			strcat(idxfn, ".bwt");
			if (!force && (fp = fopen(idxfn, "rb")) != 0) {
				fclose(fp);
				fprintf(stderr, "[E::%s] File `%s' exists. Please use `-f' to overwrite.\n", __func__, idxfn);
				return 1;
			}
		}
	}
	
	{ // read sequences
		kseq_t *seq;
		double t = cputime();
		gzFile fp;
		fp = gzopen(argv[optind], "r");
		if (fp == 0) {
			fprintf(stderr, "[E::%s] fail to open the input file.\n", __func__);
			return 1;
		}
		seq = kseq_init(fp);
		l = 0; max = 16;
		s = malloc(max);
		s[l++] = 0;
		while (kseq_read(seq) >= 0) {
			if (l + (seq->seq.l + 1) * 2 > max) {
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
		}
		kseq_destroy(seq);
		gzclose(fp);
		fprintf(stderr, "[M::%s] Load sequences in %.3f seconds.\n", __func__, cputime() - t);
	}
	
	{ // construct BWT
		int *SA = malloc(l * sizeof(int));
		double t = cputime();
		sais(s, SA, l, 6);
		fprintf(stderr, "[M::%s] Constructed the suffix array in %.3f seconds.\n", __func__, cputime() - t);
		if (check)
			fprintf(stderr, "[M::%s] SA checking returned value: %d\n", __func__, sa_check(s, SA, l));
		t = cputime();
		for (i = 0; i < l; ++i) {
			if (SA[i] == 0) SA[i] = 0;
			else SA[i] = s[SA[i] - 1];
		}
		for (i = 0; i < l; ++i) s[i] = SA[i];
		fprintf(stderr, "[M::%s] Generated BWT in %.3f seconds.\n", __func__, cputime() - t);
		free(SA);
	}

	if (!plain) {
		uint64_t len;
		int k, c;
		double t = cputime();
		rlditr_t itr;
		rld_t *e;
		e = rld_init(6, bbits);
		rld_itr_init(e, &itr, 0);
		k = 1; c = s[0];
		for (i = 1; i < l; ++i) {
			if (s[i] != c) {
				rld_enc(e, &itr, k, c);
				c = s[i];
				k = 1;
			} else ++k;
		}
		rld_enc(e, &itr, k, c);
		len = rld_enc_finish(e, &itr);
		fprintf(stderr, "[M::%s] Encoded BWT in %lld bytes in %.3f seconds\n", __func__, (unsigned long long)len, cputime() - t);
		rld_dump(e, idxfn);
		rld_destroy(e);
	} else {
		for (i = 0; i < l; ++i) putchar("$ACGTN"[s[i]]);
		putchar('\n');
	}

	free(s); free(idxfn);
	return 0;
}

int main_chkbwt(int argc, char *argv[])
{
	rld_t *e;
	rlditr_t itr;
	int i, j, l, c = 0, plain = 0;
	uint64_t *cnt, *rank, sum = 0;
	double t;
	while ((c = getopt(argc, argv, "P")) >= 0) {
		switch (c) {
			case 'P': plain = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fermi chkbwt [-P] <idxbase.bwt>\n");
		return 1;
	}
	e = rld_restore(argv[optind]);
	if (e == 0) {
		fprintf(stderr, "[E::%s] Fail to read the index file.\n", __func__);
		return 1;
	}
	cnt = alloca(e->asize * 8);
	rank = alloca(e->asize * 8);
	rld_itr_init(e, &itr, 0);
	for (i = 0; i < e->asize; ++i) cnt[i] = 0;
	t = cputime();
	while ((l = rld_dec(e, &itr, &c)) >= 0) {
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
			fprintf(stderr, "[E::%s] Different counts: %lld != %lld\n", __func__,
					(unsigned long long)cnt[j], (unsigned long long)e->mcnt[j+1]);
			exit(1);
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
	int c, n, m;
	uint64_t i, *list;
	kstring_t s;
	s.m = s.l = 0; s.s = 0;
	n = m = 0; list = 0;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
			case 'i':
				if (n == m) {
					m = m? m<<1 : 16;
					list = realloc(list, 8 * m);
				}
				list[n++] = atol(optarg); break;
				break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fermi unpack <idxbase.bwt>\n");
		return 1;
	}
	e = rld_restore(argv[optind]);
	if (n) {
		for (i = 0; (int)i < n; ++i)
			if (list[i] > 0 && list[i] < e->mcnt[1])
				print_i(e, list[i], &s);
	} else {
		for (i = 1; i < e->mcnt[1]; ++i)
			print_i(e, i, &s);
	}
	rld_destroy(e);
	free(s.s);
	return 0;
}

int main_exact(int argc, char *argv[])
{
	int c, min_match = 30;
	rld_t *e;
	kseq_t *seq;
	gzFile fp;
	while ((c = getopt(argc, argv, "m:")) >= 0) {
		switch (c) {
			case 'm': min_match = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fermi exact [-m minMatch] <idxbase.bwt> <src.fa>\n");
		return 1;
	}
	fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	e = rld_restore(argv[optind]);
	while (kseq_read(seq) >= 0) {
		seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
		printf(">%s\n", seq->name.s);
		printf("%d\n", fm6_search_forward_overlap(e, min_match, seq->seq.l, (uint8_t*)seq->seq.s));
		puts("//");
	}
	rld_destroy(e);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
