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
int main_build(int argc, char *argv[]);
int main_correct(int argc, char *argv[]);
int main_unitig(int argc, char *argv[]);
int main_clean(int argc, char *argv[]);
int main_cnt2qual(int argc, char *argv[]);
int main_seqsort(int argc, char *argv[]);
int main_pairext(int argc, char *argv[]);

int main_example(int argc, char *argv[]);

double rssmem();
double cputime();
double realtime();
void liftrlimit();

int main(int argc, char *argv[])
{
	int i, ret = 0;
	double start = realtime();
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: fermi (FERragina-Manzini Index for DNA sequences)\n");
		fprintf(stderr, "Version: %s\n", FERMI_VERSION);
		fprintf(stderr, "Contact: Heng Li <lh3@live.co.uk>\n\n");
		fprintf(stderr, "Usage:   fermi <command> [arguments]\n\n");
		fprintf(stderr, "Command: build     Generate FM-Index\n");
		fprintf(stderr, "         chkbwt    Validate the FM-Index\n");
		fprintf(stderr, "         merge     Merge multiple FM-Indexes\n");
		fprintf(stderr, "         unpack    Retrieve DNA sequences\n");
		fprintf(stderr, "         exact     Find exact matches\n");
		fprintf(stderr, "         correct   Error correction\n");
		fprintf(stderr, "         seqrank   Compute the rank of sequences\n");
		fprintf(stderr, "         unitig    Construct unitigs\n");
		fprintf(stderr, "         pairext   Unitig extension using PE reads\n");
		fprintf(stderr, "         clean     Clean the graph\n");
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
	else if (strcmp(argv[1], "cnt2qual") == 0) ret = main_cnt2qual(argc-1, argv+1);
	else if (strcmp(argv[1], "seqsort") == 0) ret = main_seqsort(argc-1, argv+1);
	else if (strcmp(argv[1], "seqrank") == 0) ret = main_seqsort(argc-1, argv+1); // an alias of seqsort
	else if (strcmp(argv[1], "unitig") == 0) ret = main_unitig(argc-1, argv+1);
	else if (strcmp(argv[1], "pairext") == 0) ret = main_pairext(argc-1, argv+1);
	else if (strcmp(argv[1], "correct") == 0) ret = main_correct(argc-1, argv+1);
	else if (strcmp(argv[1], "clean") == 0) ret = main_clean(argc-1, argv+1);
	else if (strcmp(argv[1], "splitfa") == 0) ret = main_splitfa(argc-1, argv+1);
	else if (strcmp(argv[1], "fltuniq") == 0) ret = main_fltuniq(argc-1, argv+1);
	else if (strcmp(argv[1], "trimseq") == 0) ret = main_trimseq(argc-1, argv+1);
	else if (strcmp(argv[1], "pe2cofq") == 0) ret = main_pe2cofq(argc-1, argv+1);
	else if (strcmp(argv[1], "cg2cofq") == 0) ret = main_cg2cofq(argc-1, argv+1);
	else if (strcmp(argv[1], "example") == 0) ret = main_example(argc-1, argv+1);
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
