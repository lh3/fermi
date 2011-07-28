#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <zlib.h>
#include "rld.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main_index(int argc, char *argv[])
{
	int max, l, c, bbits = 5, *SA;
	uint8_t *s;
	kseq_t *seq;
	gzFile fp;
	rldenc_t *e;

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: bwa2 index <in.fa>\n");
		return 1;
	}
	
	s = 0; l = max = 0;
	fp = gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		if (l + seq->seq.l + 1 >= max) {
			max = l + seq->seq.l + 1;
			kroundup32(max);
			s = realloc(s, max);
		}
		memcpy(s + l, seq->seq.s, seq->seq.l + 1);
		l += seq->seq.l + 1;
	}
	kseq_destroy(seq);
	gzclose(fp);

	e = rld_enc_init(6, bbits);
	free(e->cnt); free(e);
	return 0;
}
