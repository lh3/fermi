#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "fermi.h"

int main_example(int argc, char *argv[])
{
	int c, do_ec = 0, ec_k = 19, unitig_k = 30;
	char *seq, *qual;
	int64_t l;
	while ((c = getopt(argc, argv, "ek:l:")) >= 0) {
		switch (c) {
			case 'e': do_ec = 1; break;
			case 'k': ec_k = atoi(optarg); break;
			case 'l': unitig_k = atoi(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi example [-e] <in.fq>\n");
		return 1;
	}
	l = fm6_api_readseq(argv[optind], &seq, &qual);
	if (do_ec) fm6_api_correct(ec_k, l, seq, qual);
	fm6_api_writeseq(l, seq, qual);
	free(seq); free(qual);
	return 0;
}
