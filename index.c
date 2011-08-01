#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "rle6.h"
#include "rld.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt5_table[256] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

int sais(const unsigned char *T, int *SA, int n, int k);
int sa_check(const unsigned char *T, const int *SA, int n);

int main_index(int argc, char *argv[])
{
	int i, max = 0, l = 0, bbits = 5, plain = 0, use_rld = 0, check = 0;
	uint8_t *s = 0;

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
			fprintf(stderr, "Usage: bwa2 index <in.fa>\n");
			return 1;
		}
	}
	
	{ // read sequences
		kseq_t *seq;
		gzFile fp;
		fp = gzopen(argv[optind], "r");
		if (fp == 0) {
			fprintf(stderr, "[E::%s] fail to open the input file.\n", __func__);
			return 1;
		}
		seq = kseq_init(fp);
		while (kseq_read(seq) >= 0) {
			if (l + seq->seq.l + 1 >= max) {
				max = l + seq->seq.l + 1;
				kroundup32(max);
				s = realloc(s, max);
			}
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
		}
		kseq_destroy(seq);
		gzclose(fp);
		for (i = 0; i < l; ++i) s[i] = s[i] > 127? 5 : seq_nt5_table[s[i]];
	}
	
	{ // construct BWT
		int *SA = malloc(l * sizeof(int));
		sais(s, SA, l, 6);
		if (check)
			fprintf(stderr, "[M::%s] SA checking returned value: %d\n", __func__, sa_check(s, SA, l));
		for (i = 0; i < l; ++i) {
			if (SA[i] == 0) SA[i] = 0;
			else SA[i] = s[SA[i] - 1];
		}
		for (i = 0; i < l; ++i) s[i] = SA[i];
		free(SA);
	}

	if (!plain) {
		if (use_rld) {
		uint64_t len;
		int k, c;
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
		rld_index(e);
		printf("%lf\n", len/8.);
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
