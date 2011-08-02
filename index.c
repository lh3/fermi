#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "rle6.h"
#include "rld.h"
#include "exact.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int sais(const unsigned char *T, int *SA, int n, int k);
int sa_check(const unsigned char *T, const int *SA, int n);

void seq_char2nt6(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

double cputime();

int main_index(int argc, char *argv[])
{
	int i, max, l, bbits = 5, plain = 0, use_rld = 0, check = 0;
	uint8_t *s;

	{ // parse the command line
		int c;
		while ((c = getopt(argc, argv, "CPb:d")) >= 0) {
			switch (c) {
				case 'C': check = 1; break;
				case 'P': plain = 1; break;
				case 'b': bbits = atoi(optarg); break;
				case 'd': use_rld = 1; break;
			}
		}
		if (argc == optind) {
			fprintf(stderr, "Usage: fmg index <in.fa>\n");
			return 1;
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
				max = l + (seq->seq.l + 1) * 2;
				kroundup32(max);
				s = realloc(s, max);
			}
			seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
			seq_revcomp6(seq->seq.l, (uint8_t*)seq->seq.s);
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
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
		if (use_rld) {
		uint64_t len;
		int k, c;
		double t = cputime();
		rld_t *e;
		e = rld_enc_init(6, bbits);
		k = 1; c = s[0];
		for (i = 1; i < l; ++i) {
			if (s[i] != c) {
				rld_enc(e, k, c);
				c = s[i];
				k = 1;
			} else ++k;
		}
		rld_enc(e, k, c);
		len = rld_enc_finish(e);
		fprintf(stderr, "[M::%s] Encoded BWT in %lld bytes in %.3f seconds\n", __func__, len/8, cputime() - t);
		if (0) {
			int b = 3, k = 0, *SA;
			rldidx_t *r;
			r = rld_index(e);
			SA = malloc(l * sizeof(int));
			for (i = 0; i < l; ++i) {
				if (s[i] == b) ++k;
				SA[i] = k;
			}
			for (i = 0; i < l; ++i) {
				int x = rld_rank11(e, r, i, b);
				if (SA[i] != x)
					printf("fail @ %d: %d != %d\n", i, SA[i], x);
			}
			free(SA);
		}
		if (1) {
			uint64_t k, l;
			kstring_t str;
			str.l = str.m = 0; str.s = 0;
			rldidx_t *r = rld_index(e);
			printf("%lld\n", fm_backward_search(e, r, 3, (const uint8_t*)"\3\4\3", &k, &l));
			fm6_retrieve(e, r, 0, &str);
			for (i = str.l - 1; i >= 0; --i) putchar("$ACGTN"[(int)str.s[i]]); putchar('\n');
		}
		free(e->cnt); free(e);
		} else {
		uint64_t len;
		int k, c;
		rle6_t *e;
		e = rle6_enc_init();
		k = 1; c = s[0];
		for (i = 1; i < l; ++i) {
			if (s[i] != c) {
				rle6_enc(e, k, c);
				c = s[i];
				k = 1;
			} else ++k;
		}
		rle6_enc(e, k, c);
		len = rle6_enc_finish(e);
		printf("%lf\n", len/8.);
		free(e);
		}
	} else {
		for (i = 0; i < l; ++i) putchar("$ACGTN"[s[i]]);
	}

	free(s);
	return 0;
}
