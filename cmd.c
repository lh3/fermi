#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <limits.h>
#include <ctype.h>
#include "priv.h"
#include "kstring.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

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
			if (list[i] < e->mcnt[1])
				print_i(e, list[i], &s);
	} else {
		for (i = 0; i < e->mcnt[1]; ++i)
			print_i(e, i, &s);
	}
	rld_destroy(e);
	free(s.s);
	return 0;
}

static uint64_t *load_sorted(uint64_t n, const char *fn)
{
	FILE *fp;
	uint64_t *sorted;
	fp = fopen(fn, "rb");
	sorted = malloc(n * 8);
	fread(sorted, n, 8, fp);
	fclose(fp);
	return sorted;
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
		fprintf(stderr, "Usage:   fermi unitig [options] <reads.fmd>\n\n");
		fprintf(stderr, "Options: -l INT      min match [%d]\n", min_match);
		fprintf(stderr, "         -t INT      number of threads [1]\n");
		fprintf(stderr, "         -r FILE     rank file [null]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	if (fn_sorted) {
		sorted = load_sorted(e->mcnt[1], fn_sorted);
		free(fn_sorted);
	}
	fm6_unitig(e, min_match, n_threads, sorted);
	free(sorted);
	rld_destroy(e);
	return 0;
}

int main_remap(int argc, char *argv[])
{
	int c, use_mmap = 0, n_threads = 1, skip = 50, min_pcv = 0, max_dist = 1000;
	rld_t *e;
	uint64_t *sorted = 0;
	char *fn_sorted = 0;
	while ((c = getopt(argc, argv, "Ml:t:c:r:D:")) >= 0) {
		switch (c) {
			case 'l': skip = atoi(optarg); break;
			case 'M': use_mmap = 1; break;
			case 'c': min_pcv = atoi(optarg); break;
			case 't': n_threads = atoi(optarg); break;
			case 'D': max_dist = atoi(optarg); break;
			case 'r': fn_sorted = optarg; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi remap [options] <reads.fmd> <contigs.fq>\n\n");
		fprintf(stderr, "Options: -l INT      skip ending INT bases of a read pair [%d]\n", skip);
		fprintf(stderr, "         -c INT      minimum paired-end coverage [%d]\n", min_pcv);
		fprintf(stderr, "         -D INT      maximum insert size (external distance) [%d]\n", max_dist);
		fprintf(stderr, "         -r FILE     rank [null]\n");
		fprintf(stderr, "         -t INT      number of threads [1]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
	if (fn_sorted) sorted = load_sorted(e->mcnt[1], fn_sorted);
	fm6_remap(argv[optind+1], e, sorted, skip, min_pcv, max_dist, n_threads);
	free(sorted);
	rld_destroy(e);
	return 0;
}

int main_correct(int argc, char *argv[])
{
	int c, use_mmap = 0, n_threads = 1;
	rld_t *e;
	fmecopt_t opt;
	opt.w = -1; opt.min_occ = 3; opt.keep_bad = 0; opt.is_paired = 0; opt.max_corr = 0.3; opt.trim_l = 0; opt.step = 5;
	while ((c = getopt(argc, argv, "MKt:k:v:O:pC:l:s:")) >= 0) {
		switch (c) {
			case 'M': use_mmap = 1; break;
			case 'K': opt.keep_bad = 1; break;
			case 't': n_threads = atoi(optarg); break;
			case 'k': opt.w = atoi(optarg); break;
			case 'v': fm_verbose = atoi(optarg); break;
			case 'O': opt.min_occ = atoi(optarg); break;
			case 'p': opt.is_paired = 1; break;
			case 'C': opt.max_corr = atof(optarg); break;
			case 'l': opt.trim_l = atoi(optarg); break;
			case 's': opt.step = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi correct [options] <reads.fmd> <reads.fq>\n\n");
		fprintf(stderr, "Options: -k INT      k-mer length; -1 for auto [%d]\n", opt.w);
		fprintf(stderr, "         -O INT      minimum (k+1)-mer occurrences [%d]\n", opt.min_occ);
		fprintf(stderr, "         -t INT      number of threads [%d]\n", n_threads);
		fprintf(stderr, "         -C FLOAT    max fraction of corrected bases [%.2f]\n", opt.max_corr);
		fprintf(stderr, "         -l INT      trim read down to INT bp; 0 to disable [0]\n");
		fprintf(stderr, "         -s INT      step size for the jumping heuristic; 0 to disable [%d]\n", opt.step);
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
	int sbits = 3, force = 0, asize = 6, max_len = INT_MAX, no_fr = 1;
	int64_t sum_l = 0, l, max, block_size = 250000000;
	uint8_t *s;
	char *idxfn = 0;
	double t;
	rld_t *e = 0;

	{ // parse the command line
		int c;
		while ((c = getopt(argc, argv, "fb:o:i:s:l:O")) >= 0) {
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
				case 'l': max_len = atoi(optarg); break;
				case 'O': no_fr = 0; break;
			}
		}
		if (argc == optind) {
			fprintf(stderr, "\n");
			fprintf(stderr, "Usage:   fermi build [options] <in.fa>\n\n");
			fprintf(stderr, "Options: -b INT    use a small marker per 2^(INT+3) bytes [%d]\n", sbits);
			fprintf(stderr, "         -f        force to overwrite the output file (effective with -o)\n");
			fprintf(stderr, "         -i FILE   append the FM-index to the existing FILE [null]\n");
			fprintf(stderr, "         -l INT    trim read down to INT bp [inf]\n");
			fprintf(stderr, "         -o FILE   output file name [null]\n");
			fprintf(stderr, "         -O        do not trim 1bp for reads whose forward and reverse are identical\n");
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
			if (seq->seq.l > max_len)
				seq->seq.l = max_len, seq->seq.s[max_len] = 0;
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
			if (no_fr && (seq->seq.l&1) == 0) {
				int i;
				for (i = 0; i < seq->seq.l>>1; ++i)
					if (seq->seq.s[i] + seq->seq.s[seq->seq.l-1-i] != 5) break;
				if (i == seq->seq.l>>1) --seq->seq.l, seq->seq.s[seq->seq.l] = 0;
			}
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

int main_clean(int argc, char *argv[])
{
	mag_t *g;
	int c;
	magopt_t *opt;
	opt = mag_init_opt();
	while ((c = getopt(argc, argv, "ON:d:CFAl:e:i:o:R:n:w:r:S")) >= 0) {
		switch (c) {
		case 'F': opt->flag |= MOG_F_NO_AMEND; break;
		case 'C': opt->flag |= MOG_F_CLEAN; break;
		case 'A': opt->flag |= MOG_F_AGGRESSIVE; break;
		case 'O': opt->flag |= MOG_F_READ_ORI; break;
		case 'S': opt->flag |= MOG_F_NO_SIMPL; break;
		case 'd': opt->min_dratio0 = atof(optarg); break;
		case 'N': opt->max_arc  = atoi(optarg); break;
		case 'l': opt->min_elen = atoi(optarg); break;
		case 'e': opt->min_ensr = atoi(optarg); break;
		case 'i': opt->min_insr = atoi(optarg); break;
		case 'o': opt->min_ovlp = atoi(optarg); break;
		case 'n': opt->n_iter   = atoi(optarg); break;
		case 'R': opt->min_dratio1 = atof(optarg); break;
		case 'w': opt->max_bcov = atof(optarg); break;
		case 'r': opt->max_bfrac= atof(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi clean [options] <in.mog>\n\n");
		fprintf(stderr, "Options: -N INT      read maximum INT neighbors per node [%d]\n", opt->max_arc);
		fprintf(stderr, "         -d FLOAT    drop a neighbor if relative overlap ratio below FLOAT [%.2f]\n\n", opt->min_dratio0); 
		fprintf(stderr, "         -C          clean the graph\n");
		fprintf(stderr, "         -l INT      minimum tip length [%d]\n", opt->min_elen);
		fprintf(stderr, "         -e INT      minimum tip read count [%d]\n", opt->min_ensr);
		fprintf(stderr, "         -i INT      minimum internal unitig read count [%d]\n", opt->min_insr);
		fprintf(stderr, "         -o INT      minimum overlap [%d]\n", opt->min_ovlp);
		fprintf(stderr, "         -R FLOAT    minimum relative overlap ratio [%.2f]\n", opt->min_dratio1);
		fprintf(stderr, "         -n INT      number of iterations [%d]\n", opt->n_iter);
		fprintf(stderr, "         -A          aggressive bubble popping\n");
		fprintf(stderr, "         -S          skip bubble simplification\n");
		fprintf(stderr, "         -w FLOAT    minimum coverage to keep a bubble [%.2f]\n", opt->max_bcov);
		fprintf(stderr, "         -r FLOAT    minimum fraction to keep a bubble [%.2f]\n", opt->max_bfrac);
		fprintf(stderr, "\n");
		return 1;
	}
	g = mag_g_read(argv[optind], opt);
	mag_g_clean(g, opt);
	mag_g_print(g);
	mag_g_destroy(g);
	free(opt);
	return 0;
}

int main_scaf(int argc, char *argv[])
{
	int c, n_threads = 1;
	rld_t *e;
	fmscafopt_t opt;
	opt.min_supp = 5; opt.pr_links = 0; opt.a_thres = 20.; opt.p_thres = 1e-20;
	while ((c = getopt(argc, argv, "m:t:Pea:p:")) >= 0) {
		switch (c) {
			case 't': n_threads = atoi(optarg); break;
			case 'm': opt.min_supp = atoi(optarg); break;
			case 'P': opt.pr_links = 1; break;
			case 'a': opt.a_thres = atof(optarg); break;
			case 'p': opt.p_thres = atof(optarg); break;
		}
	}
	if (optind + 4 > argc) {
		fprintf(stderr, "\nUsage:   fermi scaf <in.fmd> <in.remapped.mag> <avg> <std>\n\n");
		fprintf(stderr, "Options: -t INT     number of threads [1]\n");
		fprintf(stderr, "         -m INT     minimum number of supporting reads [%d]\n", opt.min_supp);
		fprintf(stderr, "         -P         print the links between unitigs\n\n");
		return 1;
	}
	opt.avg = atof(argv[optind+2]); opt.std = atof(argv[optind+3]);
	e = rld_restore(argv[optind]);
	mag_scaf_core(e, argv[optind+1], &opt, n_threads);
	rld_destroy(e);
	return 0;
}

int main_contrast(int argc, char *argv[])
{
	extern void fm6_contrast(rld_t *const e[2], int k, int min_occ, int n_threads, uint64_t *sub[2]);
	extern int64_t fm6_sub_conv(int64_t n_seqs, uint64_t *sub, const uint64_t *rank);

	int c, min_occ = 3, k = 55, n_threads = 1, i;
	uint64_t *sub[2], n_seqs[2], *rank;
	rld_t *e[2];

	while ((c = getopt(argc, argv, "k:o:t:")) >= 0) {
		switch (c) {
		case 'k': k = atoi(optarg); break;
		case 'o': min_occ = atoi(optarg); break;
		case 't': n_threads = atoi(optarg); break;
		}
	}
	if (optind + 6 > argc) {
		fprintf(stderr, "\nUsage:   fermi contrast <idx1.fmd> <idx1.rank> <1-2.sub> <idx2.fmd> <idx2.rank> <2-1.sub>\n\n");
		fprintf(stderr, "Options: -o INT    minimum occurrence [%d]\n", min_occ);
		fprintf(stderr, "         -t INT    number of threads [%d]\n", n_threads);
		fprintf(stderr, "         -k INT    k-mer length [%d]\n\n", k);
		return 1;
	}

	e[0] = rld_restore(argv[optind+0]);
	e[1] = rld_restore(argv[optind+3]);
	fm6_contrast(e, k, min_occ, n_threads, sub);
	n_seqs[0] = e[0]->mcnt[1];
	n_seqs[1] = e[1]->mcnt[1];
	rld_destroy(e[0]);
	rld_destroy(e[1]);

	rank = malloc((n_seqs[0] > n_seqs[1]? n_seqs[0] : n_seqs[1]) * 8);
	for (i = 0; i < 2; ++i) {
		FILE *fp;
		uint64_t n_sel;
		fp = fopen(argv[optind+i*3+1], "rb");
		fread(rank, 8, n_seqs[i], fp);
		fclose(fp);
		n_sel = fm6_sub_conv(n_seqs[i], sub[i], rank);
		fprintf(stderr, "[M::%s] %ld reads selected from %s\n", __func__, (long)n_sel, argv[optind+i*3]);
		fp = fopen(argv[optind+i*3+2], "wb");
		fwrite(&n_seqs[i], 8, 1, fp);
		fwrite(sub[i], 8, (n_seqs[i] + 63) / 64, fp);
		fclose(fp);
		free(sub[i]);
	}
	free(rank);
	return 0;
}

int main_sub(int argc, char *argv[])
{
	int c, n_threads = 1, is_comp = 0;
	rld_t *e;
	uint64_t n_seqs, *sub;
	FILE *fp;
	while ((c = getopt(argc, argv, "ct:")) >= 0) {
		switch (c) {
		case 't': n_threads = atoi(optarg); break;
		case 'c': is_comp = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fermi sub [-c] [-t nThreads] <in.fmd> <array.bits>\n");
		return 1;
	}
	e = rld_restore(argv[optind]);
	fp = fopen(argv[optind+1], "rb");
	fread(&n_seqs, 8, 1, fp);
	if (n_seqs != e->mcnt[1]) {
		fprintf(stderr, "[E::%s] unmatched index and the bit array\n", __func__);
		rld_destroy(e);
		return 1;
	}
	sub = malloc((n_seqs + 63) / 64 * 8);
	fread(sub, 8, (n_seqs + 63) / 64, fp);
	fclose(fp);
	e = fm_sub(e, sub, n_threads, is_comp);
	free(sub);
	rld_dump(e, "-");
	rld_destroy(e);
	return 0;
}

int main_recode(int argc, char *argv[])
{
	rld_t *e;
	if (argc == 1) {
		fprintf(stderr, "Usage: fermi recode <in.rld>\n");
		return 1;
	}
	e = rld_restore(argv[1]);
	rld_dump(e, "-");
	rld_destroy(e);
	return 0;
}

static inline uint64_t popcount64(uint64_t y)
{
	y = y - ((y >> 1) & 0x5555555555555555ull);
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static uint64_t count_bits(uint64_t len, uint64_t *x)
{
	uint64_t k, cnt;
	for (k = cnt = 0; k < len; ++k)
		cnt += popcount64(x[k]);
	return cnt;
}

static uint64_t *read_sub(const char *fn, uint64_t *len)
{
	uint64_t *sub, end;
	FILE *fp;
	*len = 0;
	if ((fp = fopen(fn, "rb")) == 0) return 0;
	fread(len, 8, 1, fp);
	end = (*len + 63) / 64;
	sub = calloc(end, 8);
	fread(sub, 8, end, fp);
	fclose(fp);
	fprintf(stderr, "[M::%s] loaded file `%s' containing %ld bits\n", __func__, fn, (long)count_bits(end, sub));
	return sub;
}

int main_bitand(int argc, char *argv[])
{
	int i;
	uint64_t len0, *sub0, k, end;
	if (argc < 3) {
		fprintf(stderr, "Usage: fermi bitand <in1.bit> <in2.bit> [...]\n");
		return 1;
	}
	sub0 = read_sub(argv[1], &len0);
	end = (len0 + 63) / 64;
	for (i = 2; i < argc; ++i) {
		uint64_t len1, *sub1;
		sub1 = read_sub(argv[i], &len1);
		if (len1 != len0) {
			fprintf(stderr, "[E::%s] unequal array length\n", __func__);
			return 1; // memory leak
		}
		for (k = 0; k < end; ++k)
			sub0[k] &= sub1[k];
		free(sub1);
	}
	fprintf(stderr, "[M::%s] the output contains %ld bits\n", __func__, (long)count_bits(end, sub0));
	fwrite(&len0, 8, 1, stdout);
	fwrite(sub0, 8, end, stdout);
	free(sub0);
	return 0;
}
