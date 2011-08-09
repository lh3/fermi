#include <string.h>
#include <stdio.h>
#include "fermi.h"

int main_index(int argc, char *argv[]);
int main_chkbwt(int argc, char *argv[]);
int main_unpack(int argc, char *argv[]);
int main_exact(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: fermi (FERragina-Manzini Index for DNA sequences)\n");
		fprintf(stderr, "Version: %s\n", FERMI_VERSION);
		fprintf(stderr, "Contact: Heng Li <lh3@live.co.uk>\n\n");
		fprintf(stderr, "Usage:   fermi <command> [arguments]\n\n");
		fprintf(stderr, "Command: index     Generated FM-Index for sequences shorter than 1G\n");
		fprintf(stderr, "         chkbwt    Validate the FM-Index\n");
		fprintf(stderr, "         unpack    Retrieve DNA sequences\n");
		fprintf(stderr, "         exact     Find exact matches\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "index") == 0) return main_index(argc-1, argv+1);
	else if (strcmp(argv[1], "chkbwt") == 0) return main_chkbwt(argc-1, argv+1);
	else if (strcmp(argv[1], "unpack") == 0) return main_unpack(argc-1, argv+1);
	else if (strcmp(argv[1], "exact") == 0) return main_exact(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unrecognized command.\n", __func__);
		return -1;
	}
	return 0;
}
