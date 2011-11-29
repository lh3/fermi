#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include "utils.h"
#include "kstring.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void seq_char2nt6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
}

void seq_reverse(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		s[l-1-i] = s[i]; s[i] = tmp;
	}
}

void seq_comp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i)
		s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

void seq_revcomp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

/**********************
 * Command-line tools *
 **********************/

int main_splitfa(int argc, char *argv[])
{
	int64_t n_seqs = 0;
	int i, n_files = 8;
	gzFile fp, *out;
	kseq_t *seq;
	char *str;
	kstring_t *ss;

	if (argc < 3) {
		fprintf(stderr, "Usage: fermi splitfa <in.fq> <out.prefix> [%d]\n", n_files);
		return 1;
	}
	if (argc >= 4) n_files = atoi(argv[3]);
	out = calloc(n_files, sizeof(gzFile));
	str = calloc(strlen(argv[2]) + 20, 1);
	ss = calloc(n_files, sizeof(kstring_t));
	for (i = 0; i < n_files; ++i) {
		sprintf(str, "%s.%.4d.fq.gz", argv[2], i);
		out[i] = gzopen(str, "wb1");
	}
	fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		i = (n_seqs>>1) % n_files;
		kputc(seq->qual.l? '@' : '>', &ss[i]); kputsn(seq->name.s, seq->name.l, &ss[i]); kputc('\n', &ss[i]);
		kputsn(seq->seq.s, seq->seq.l, &ss[i]); kputc('\n', &ss[i]);
		if (seq->qual.l) {
			kputsn("+\n", 2, &ss[i]); kputsn(seq->qual.s, seq->qual.l, &ss[i]); kputc('\n', &ss[i]);
		}
		if (ss[i].l > 64000) {
			gzwrite(out[i], ss[i].s, ss[i].l);
			ss[i].l = 0;
		}
		++n_seqs;
	}
	for (i = 0; i < n_files; ++i) {
		gzwrite(out[i], ss[i].s, ss[i].l);
		gzclose(out[i]);
		free(ss[i].s);
	}
	free(out); free(ss); free(str);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int main_fltuniq(int argc, char *argv[])
{
	int c, k = 17;
	uint64_t *flags, mask;
	kstring_t out;
	gzFile fp;
	kseq_t *seq;

	while ((c = getopt(argc, argv, "k:")) >= 0) {
		switch (c) {
			case 'k': k = atoi(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi fltuniq [-k kmer=17] <in.fa>\n");
		return 1;
	}

	out.l = out.m = 0; out.s = 0;
	fp = gzopen(argv[optind], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, argv[optind]);
		return 1;
	}
	seq = kseq_init(fp);
	flags = xcalloc(1ULL<<k*2>>5, 8);
	mask = (1ULL<<k*2) - 1;
	fprintf(stderr, "[M::%s] building the hash table...\n", __func__);
	while (kseq_read(seq) >= 0) {
		int i, l;
		uint64_t z;
		for (i = l = 0, z = 0; i < seq->seq.l; ++i) {
			c = seq_nt6_table[(int)seq->seq.s[i]] - 1;
			if (c < 4) {
				++l, z = (z<<2 | c) & mask;
				if (l >= k)
					flags[z>>5] |= ((flags[z>>5]>>((z&31)<<1)&3) == 0? 1ULL : 3ULL) << ((z&31)<<1);
			} else l = 0, z = 0;
		}
	}
	kseq_destroy(seq);
	gzrewind(fp);
	seq = kseq_init(fp);
	fprintf(stderr, "[M::%s] filtering the reads...\n", __func__);
	while (kseq_read(seq) >= 0) {
		int i, l;
		uint64_t z;
		for (i = l = 0, z = 0; i < seq->seq.l; ++i) {
			c = seq_nt6_table[(int)seq->seq.s[i]] - 1;
			if (c < 4) {
				++l, z = (z<<2 | c) & mask;
				if (l >= k && (flags[z>>5]>>((z&31)<<1)&3) != 3) break;
			} else break;
		}
		if (i == seq->seq.l) { // unfiltered
			out.l = 0;
			kputc('@', &out); kputsn(seq->name.s, seq->name.l, &out); kputc('\n', &out);
			kputsn(seq->seq.s, seq->seq.l, &out); kputc('\n', &out);
			if (seq->qual.l) {
				kputsn("+\n", 2, &out);
				kputsn(seq->qual.s, seq->qual.l, &out); kputc('\n', &out);
			}
			fputs(out.s, stdout);
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
	free(flags); free(out.s);
	return 0;
}
