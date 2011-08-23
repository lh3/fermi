#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include "fermi.h"

int main_chkbwt(int argc, char *argv[]);
int main_unpack(int argc, char *argv[]);
int main_exact(int argc, char *argv[]);
int main_merge(int argc, char *argv[]);
int main_build(int argc, char *argv[]);
int main_strlen(int argc, char *argv[]);

double rssmem();
double cputime();
double realtime();
void liftrlimit();

int main(int argc, char *argv[])
{
	int ret = 0;
	double start = realtime();
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: fermi (FERragina-Manzini Index for DNA sequences)\n");
		fprintf(stderr, "Version: %s\n", FERMI_VERSION);
		fprintf(stderr, "Contact: Heng Li <lh3@live.co.uk>\n\n");
		fprintf(stderr, "Usage:   fermi <command> [arguments]\n\n");
		fprintf(stderr, "Command: strlen    Total number of symbols to index\n");
		fprintf(stderr, "         build     Generate FM-Index\n");
		fprintf(stderr, "         chkbwt    Validate the FM-Index\n");
		fprintf(stderr, "         merge     Merge two FM-Indexes\n");
		fprintf(stderr, "         unpack    Retrieve DNA sequences\n");
		fprintf(stderr, "         exact     Find exact matches\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "build") == 0) ret = main_build(argc-1, argv+1);
	else if (strcmp(argv[1], "chkbwt") == 0) ret = main_chkbwt(argc-1, argv+1);
	else if (strcmp(argv[1], "unpack") == 0) ret = main_unpack(argc-1, argv+1);
	else if (strcmp(argv[1], "exact") == 0) ret = main_exact(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) ret = main_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "strlen") == 0) ret = main_strlen(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unrecognized command.\n", __func__);
		return -1;
	}
	if (ret == 0 && fm_verbose >= 3)
		fprintf(stderr, "[M::%s] Real time: %.3f sec; CPU: %.3f sec; RSS: %.3f MB\n", __func__, realtime() - start, cputime(), rssmem());
	return ret;
}
