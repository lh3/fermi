#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include "fermi.h"

int main_splitfa(int argc, char *argv[]);
int main_fltuniq(int argc, char *argv[]);
int main_cg2cofq(int argc, char *argv[]);
int main_pe2cofq(int argc, char *argv[]);
int main_trimseq(int argc, char *argv[]);

int main_chkbwt(int argc, char *argv[]);
int main_unpack(int argc, char *argv[]);
int main_exact(int argc, char *argv[]);
int main_merge(int argc, char *argv[]);
int main_sub(int argc, char *argv[]);
int main_build(int argc, char *argv[]);
int main_correct(int argc, char *argv[]);
int main_unitig(int argc, char *argv[]);
int main_clean(int argc, char *argv[]);
int main_cnt2qual(int argc, char *argv[]);
int main_seqsort(int argc, char *argv[]);
int main_remap(int argc, char *argv[]);
int main_scaf(int argc, char *argv[]);
int main_contrast(int argc, char *argv[]);
int main_bitand(int argc, char *argv[]);
int main_recode(int argc, char *argv[]);

int main_ropebwt(int argc, char *argv[]);
int main_example(int argc, char *argv[]);

double rssmem();
double cputime();
double realtime();
void liftrlimit();
/*
#include "rld.h"
int main_test(int argc, char *argv[])
{
	extern uint64_t fm_multi_backward_search(int n, rld_t *const*e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end);
	int i, l;
	uint64_t beg, end;
	l = strlen(argv[1]);
	for (i = 0; i < l; ++i)
		argv[1][i] = seq_nt6_table[(int)argv[1][i]];
	if (argc > 3) {
		int n = argc - 2;
		rld_t **e;
		e = alloca(n * sizeof(void*));
		for (i = 2; i < argc; ++i)
			e[i - 2] = rld_restore(argv[i]);
		fm_multi_backward_search(n, e, l, (uint8_t*)argv[1], &beg, &end);
		printf("%lld,%lld\n", beg, end);
	} else {
		rld_t *e;
		e = rld_restore(argv[2]);
		fm_backward_search(e, l, (uint8_t*)argv[1], &beg, &end);
		printf("%lld,%lld\n", beg, end);
	}
	return 0;
}
*/
int main(int argc, char *argv[])
{
	int i, ret = 0;
	double start = realtime();
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: fermi (FERragina-Manzini Index for DNA sequences)\n");
		fprintf(stderr, "Version: %s\n", FERMI_VERSION);
		fprintf(stderr, "Contact: <http://hengli.uservoice.com>\n\n");
		fprintf(stderr, "Usage:   fermi <command> [arguments]\n\n");
		fprintf(stderr, "Command: build     Generate FM-Index\n");
		fprintf(stderr, "         ropebwt   Alternative algorithms for constructing FM-index\n");
		fprintf(stderr, "         chkbwt    Validate the FM-Index\n");
		fprintf(stderr, "         merge     Merge multiple FM-Indices\n");
		fprintf(stderr, "         unpack    Retrieve DNA sequences\n");
		fprintf(stderr, "         exact     Find exact matches\n");
		fprintf(stderr, "         correct   Error correction\n");
		fprintf(stderr, "         seqrank   Compute the rank of sequences\n");
		fprintf(stderr, "         unitig    Construct unitigs\n");
		fprintf(stderr, "         clean     Clean the graph\n");
		fprintf(stderr, "         remap     Compute the coverage and PE coverage\n");
		fprintf(stderr, "         scaf      Generate scaftigs\n");
		fprintf(stderr, "         contrast  Compare two FMD-indices\n");
		fprintf(stderr, "         bitand    Compute the intersection of bit arrays\n");
		fprintf(stderr, "         sub       Extract sub-index with a bit array\n");
		fprintf(stderr, "         recode    Recode FM-Index\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "         splitfa   Split a FASTA/Q file\n");
		fprintf(stderr, "         trimseq   Trim a FASTA/Q file\n");
		fprintf(stderr, "         fltuniq   Filter out reads containing unique mer\n");
		fprintf(stderr, "         pe2cofq   Convert split pefq to collate fastq\n");
		fprintf(stderr, "         cg2cofq   Convert cgfq to collate fastq (deprecated)\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "         example   Light-weight assembly using fermi high-level APIs\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "build") == 0) ret = main_build(argc-1, argv+1);
	else if (strcmp(argv[1], "chkbwt") == 0) ret = main_chkbwt(argc-1, argv+1);
	else if (strcmp(argv[1], "unpack") == 0) ret = main_unpack(argc-1, argv+1);
	else if (strcmp(argv[1], "exact") == 0) ret = main_exact(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) ret = main_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "sub") == 0) ret = main_sub(argc-1, argv+1);
	else if (strcmp(argv[1], "cnt2qual") == 0) ret = main_cnt2qual(argc-1, argv+1);
	else if (strcmp(argv[1], "seqsort") == 0) ret = main_seqsort(argc-1, argv+1);
	else if (strcmp(argv[1], "seqrank") == 0) ret = main_seqsort(argc-1, argv+1); // an alias of seqsort
	else if (strcmp(argv[1], "unitig") == 0) ret = main_unitig(argc-1, argv+1);
	else if (strcmp(argv[1], "remap") == 0) ret = main_remap(argc-1, argv+1);
	else if (strcmp(argv[1], "scaf") == 0) ret = main_scaf(argc-1, argv+1);
	else if (strcmp(argv[1], "correct") == 0) ret = main_correct(argc-1, argv+1);
	else if (strcmp(argv[1], "clean") == 0) ret = main_clean(argc-1, argv+1);
	else if (strcmp(argv[1], "recode") == 0) ret = main_recode(argc-1, argv+1);
	else if (strcmp(argv[1], "splitfa") == 0) ret = main_splitfa(argc-1, argv+1);
	else if (strcmp(argv[1], "fltuniq") == 0) ret = main_fltuniq(argc-1, argv+1);
	else if (strcmp(argv[1], "trimseq") == 0) ret = main_trimseq(argc-1, argv+1);
	else if (strcmp(argv[1], "pe2cofq") == 0) ret = main_pe2cofq(argc-1, argv+1);
	else if (strcmp(argv[1], "cg2cofq") == 0) ret = main_cg2cofq(argc-1, argv+1);
	else if (strcmp(argv[1], "example") == 0) ret = main_example(argc-1, argv+1);
	else if (strcmp(argv[1], "contrast") == 0) ret = main_contrast(argc-1, argv+1);
	else if (strcmp(argv[1], "bitand") == 0) ret = main_bitand(argc-1, argv+1);
	else if (strcmp(argv[1], "ropebwt") == 0) ret = main_ropebwt(argc-1, argv+1);
//	else if (strcmp(argv[1], "test") == 0) ret = main_test(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unrecognized command.\n", __func__);
		return -1;
	}
	if (ret == 0 && fm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, FERMI_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; RSS: %.3f MB\n", __func__, realtime() - start, cputime(), rssmem());
	}
	return ret;
}
