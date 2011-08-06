#include <string.h>
#include <stdio.h>

int main_index(int argc, char *argv[]);
int main_chkbwt(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "!!!\n");
		return 1;
	}
	if (strcmp(argv[1], "index") == 0) return main_index(argc-1, argv+1);
	else if (strcmp(argv[1], "chkbwt") == 0) return main_chkbwt(argc-1, argv+1);
	else {
		fprintf(stderr, "[E::%s] unrecognized command.\n", __func__);
		return -1;
	}
	return 0;
}
