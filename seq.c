#include <zlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include "utils.h"
#include "kstring.h"
#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

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

static void write_seq(const kseq_t *seq, kstring_t *out)
{
	kputc(seq->qual.l? '@' : '>', out);
	kputsn(seq->name.s, seq->name.l, out);
	if (seq->comment.l) {
		kputc(' ', out);
		kputsn(seq->comment.s, seq->comment.l, out);
	}
	kputc('\n', out);
	kputsn(seq->seq.s, seq->seq.l, out);
	if (seq->qual.l) {
		kputsn("\n+\n", 3, out);
		kputsn(seq->qual.s, seq->qual.l, out);
	}
	kputc('\n', out);
}

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
		write_seq(seq, &ss[i]);
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
	int c, k = 0;
	uint64_t *flags, mask;
	kstring_t out, prev_name;
	gzFile fp;
	kseq_t *seq;

	while ((c = getopt(argc, argv, "k:")) >= 0) {
		switch (c) {
			case 'k': k = atoi(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi fltuniq <in.fa>\n");
		return 1;
	}

	if (k == 0) { // compute the k-mer length based on the input file size
		FILE *fp;
		if ((fp = fopen(argv[optind], "rb")) == 0) {
			fprintf(stderr, "[E::%s] fail to open the input file\n", __func__);
			return 1;
		} else {
			fseek(fp, 0, SEEK_END);
			k = (int)(log(ftell(fp)) / log(4) + 1.499);
			fclose(fp);
			if (k > 18) k = 18;
			if (k < 15) k = 15;
			fprintf(stderr, "[M::%s] set the k-mer size as %d\n", __func__, k);
		}
	}

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
	out.l = out.m = 0; out.s = 0;
	prev_name.l = prev_name.m = 0; prev_name.s = 0;
	fprintf(stderr, "[M::%s] filtering the reads...\n", __func__);
	while (kseq_read(seq) >= 0) {
		int i, l, is_paired;
		uint64_t z;
		is_paired = (prev_name.l && strcmp(prev_name.s, seq->name.s) == 0);
		if (is_paired) { // see trimseq() for more comments
			if (out.l == 0) continue;
		} else {
			if (out.l) fputs(out.s, stdout);
			out.l = 0;
		}
		for (i = l = 0, z = 0; i < seq->seq.l; ++i) {
			c = seq_nt6_table[(int)seq->seq.s[i]] - 1;
			if (c < 4) {
				++l, z = (z<<2 | c) & mask;
				if (l >= k && (flags[z>>5]>>((z&31)<<1)&3) != 3) break;
			} else break;
		}
		if (i == seq->seq.l) write_seq(seq, &out);
		else if (is_paired) out.l = 0;
		prev_name.l = 0;
		kputsn(seq->name.s, seq->name.l, &prev_name);
	}
	if (out.l) fputs(out.s, stdout);

	kseq_destroy(seq);
	gzclose(fp);
	free(flags); free(out.s);
	return 0;
}

int main_cg2cofq(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	kstring_t str;

	if (argc == 1) {
		fprintf(stderr, "Usage: fermi cg2cofq <in.cgfq>\n");
		return 1;
	}
	str.l = str.m = 0; str.s = 0;
	fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i;
		str.l = 0;
		ks_resize(&str, seq->name.l + seq->seq.l + seq->qual.l + 32);
		kputc(seq->qual.l? '@' : '>', &str); kputsn(seq->name.s, seq->name.l, &str); kputc('\n', &str);
		for (i = 0; i < seq->seq.l; ++i) 
			if (!isalpha(seq->seq.s[i])) break;
		kputsn(seq->seq.s, i, &str);
		if (seq->qual.l) {
			kputsn("\n+\n", 3, &str);
			kputsn(seq->qual.s, i, &str);
		}
		kputc('\n', &str);
		for (; seq->seq.l; ++i)
			if (isalpha(seq->seq.s[i])) break;
		if (i != seq->seq.l) {
			kputc(seq->qual.l? '@' : '>', &str); kputsn(seq->name.s, seq->name.l, &str); kputc('\n', &str);
			kputsn(seq->seq.s + i, seq->seq.l - i, &str);
			if (seq->qual.l) {
				kputsn("\n+\n", 3, &str);
				kputsn(seq->qual.s + i, seq->qual.l - i, &str);
			}
			kputc('\n', &str);
		}
		fputs(str.s, stdout);
	}
	kseq_destroy(seq);
	gzclose(fp);
	free(str.s);
	return 0;
}

int main_pe2cofq(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *seq[2];
	kstring_t str;

	if (argc < 3) {
		fprintf(stderr, "Usage: fermi pe2cofq <in1.fq> <in2.fq>\n");
		return 1;
	}
	str.l = str.m = 0; str.s = 0;
	fp1 = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	fp2 = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
	seq[0] = kseq_init(fp1);
	seq[1] = kseq_init(fp2);
	while (kseq_read(seq[0]) >= 0) {
		if (kseq_read(seq[1]) < 0) break; // one file ends
		str.l = 0;
		if (seq[0]->name.l > 2 && seq[0]->name.s[seq[0]->name.l-2] == '/' && isdigit(seq[0]->name.s[seq[0]->name.l-1]))
			seq[0]->name.s[(seq[0]->name.l -= 2)] = 0; // trim tailing "/[0-9]$"
		seq[1]->name.l = 0;
		kputsn(seq[0]->name.s, seq[0]->name.l, &seq[1]->name); // make sure two ends having the same name
		write_seq(seq[0], &str);
		write_seq(seq[1], &str);
		fputs(str.s, stdout);
	}
	kseq_destroy(seq[0]); gzclose(fp1);
	kseq_destroy(seq[1]); gzclose(fp2);
	free(str.s);
	return 0;
}

int main_trimseq(int argc, char *argv[])
{
	int c, min_l = 20, min_q = 3, drop_ambi = 1;
	gzFile fp;
	kseq_t *seq;
	kstring_t prev_name, str;

	while ((c = getopt(argc, argv, "q:Nl:")) >= 0) {
		switch (c) {
			case 'q': min_q = atoi(optarg); break;
			case 'l': min_l = atoi(optarg); break;
			case 'N': drop_ambi = 0; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fermi trimseq [-N] [-q qual=%d] [-l minLen=%d] <in.fq>\n", min_q, min_l);
		return 1;
	}
	str.l = str.m = 0; str.s = 0;
	prev_name.l = prev_name.m = 0; prev_name.s = 0;
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i, is_paired = 0, left, right, drop = 0;
		if (seq->name.l == prev_name.l && prev_name.l) { // test pairing
			if (strncmp(seq->name.s, prev_name.s, seq->name.l - 1) == 0) {
				int c2 = seq->name.s[prev_name.l-1], c1 = prev_name.s[prev_name.l-1];
				if (c1 == c2) is_paired = 1;
				else if (prev_name.l >= 2 && prev_name.s[prev_name.l-2] == '/') {
					if (isdigit(c1) && isdigit(c2))
						is_paired = 1;
				}
			}
		}
		if (is_paired) {
			if (str.l == 0) continue; // do not process this sequence
		} else { // output the previous sequence(s)
			if (str.l) fputs(str.s, stdout);
			str.l = 0;
		}
		// process the current sequence
		left = 0; right = seq->seq.l;
		if (min_q > 0 && seq->qual.l) { // trim
			int s, max, max_i;
			// trim from 3'-end
			for (i = right - 1, max = s = 0, max_i = right; i >= left; --i) {
				s += min_q - (seq->qual.s[i] - 33);
				if (s < 0) break;
				if (max < s) max = s, max_i = i;
			}
			right = max_i;
			// trim from 5'-end
			for (i = 0, max = s = 0, max_i = -1; i < right; ++i) {
				s += min_q - (seq->qual.s[i] - 33);
				if (s < 0) break;
				if (max < s) max = s, max_i = i;
			}
			left = max_i + 1;
			if (right - left < min_l) drop = 1;
		}
		if (!drop && drop_ambi) {
			for (i = left; i < right; ++i)
				if (seq_nt6_table[(int)seq->seq.s[i]] >= 5)
					break;
			if (i != right) drop = 1;
		}
		if (!drop) {
			if (left) {
				memmove(seq->seq.s, seq->seq.s + left, right - left);
				if (seq->qual.l)
					memmove(seq->qual.s, seq->qual.s + left, right - left);
			}
			seq->seq.l = right - left;
			if (seq->qual.l) seq->qual.l = right - left;
			write_seq(seq, &str);
		} else if (is_paired) str.l = 0;
		prev_name.l = 0;
		kputsn(seq->name.s, seq->name.l, &prev_name);
	}
	if (str.l) fputs(str.s, stdout);
	kseq_destroy(seq);
	gzclose(fp);
	free(str.s); free(prev_name.s);
	return 0;
}

/*******************
 * High-level APIs *
 *******************/

#include "kvec.h"
#include "ksort.h"
KSORT_INIT_GENERIC(int)

#define DEFAULT_QUAL 20

int64_t fm6_api_readseq(const char *fn, char **_seq, char **_qual)
{
	gzFile fp;
	kseq_t *kseq;
	kstring_t seq, qual;
	seq.l = seq.m = qual.l = qual.m = 0; seq.s = qual.s = 0;
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	kseq = kseq_init(fp);
	while (kseq_read(kseq) >= 0) {
		kputsn(kseq->seq.s, kseq->seq.l + 1, &seq);
		if (kseq->qual.l == 0) {
			int i;
			ks_resize(&qual, seq.m);
			for (i = 0; i < kseq->seq.l; ++i)
				kputc(DEFAULT_QUAL + 33, &qual);
			kputc(0, &qual);
		} else kputsn(kseq->qual.s, kseq->qual.l + 1, &qual);
	}
	assert(seq.l == qual.l);
	kseq_destroy(kseq);
	gzclose(fp);
	*_seq = seq.s; *_qual = qual.s;
	return seq.l;
}

void fm6_api_writeseq(int64_t l, char *seq, char *qual)
{
	kstring_t s, q;
	int64_t i;
	s.l = s.m = q.l = q.m = 0; s.s = q.s = 0;
	for (i = 0; i < l; ++i) {
		int c = seq[i];
		if (c == 0) {
			printf("@%ld\n", (long)i);
			puts(s.s); puts("+"); puts(q.s);
			s.l = q.l = 0;
		} else {
			if (c < 6) c = "$ACGTN"[c];
			kputc(c, &s);
			kputc(qual[i], &q);
		}
	}
	free(s.s); free(q.s);
}

int fm6_api_seqlen(int64_t l, const char *seq, double quantile)
{
	int64_t beg, i, cnt, j;
	int *len, ret;
	for (i = cnt = 0; i < l; ++i)
		if (seq[i] == 0) ++cnt;
	len = malloc(cnt * sizeof(int));
	for (beg = i = j = 0; i < l; ++i)
		if (seq[i] == 0)
			len[j++] = i - beg, beg = i + 1;
	assert(j == cnt);
	ret = ks_ksmall(int, cnt, len, cnt * quantile);
	free(len);
	return ret;
}
