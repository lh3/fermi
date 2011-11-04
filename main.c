#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include "fermi.h"

int main_chkbwt(int argc, char *argv[]);
int main_unpack(int argc, char *argv[]);
int main_exact(int argc, char *argv[]);
int main_merge(int argc, char *argv[]);
int main_build(int argc, char *argv[]);
int main_splitfa(int argc, char *argv[]);
int main_correct(int argc, char *argv[]);
int main_unitig(int argc, char *argv[]);
int main_clean(int argc, char *argv[]);
int main_cnt2qual(int argc, char *argv[]);

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
		fprintf(stderr, "         unitig    Unitig\n");
		fprintf(stderr, "         clean     Clean the graph\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "build") == 0) ret = main_build(argc-1, argv+1);
	else if (strcmp(argv[1], "chkbwt") == 0) ret = main_chkbwt(argc-1, argv+1);
	else if (strcmp(argv[1], "unpack") == 0) ret = main_unpack(argc-1, argv+1);
	else if (strcmp(argv[1], "exact") == 0) ret = main_exact(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) ret = main_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "splitfa") == 0) ret = main_splitfa(argc-1, argv+1);
	else if (strcmp(argv[1], "cnt2qual") == 0) ret = main_cnt2qual(argc-1, argv+1);
	else if (strcmp(argv[1], "unitig") == 0) ret = main_unitig(argc-1, argv+1);
	else if (strcmp(argv[1], "correct") == 0) ret = main_correct(argc-1, argv+1);
	else if (strcmp(argv[1], "clean") == 0) ret = main_clean(argc-1, argv+1);
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
