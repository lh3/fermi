#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bprope6.h"
#include "bcr.h"
#include "rld.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char seq_nt6_table[128];

enum algo_e { BPR, RBR, BCR };

#define FLAG_FOR 0x1
#define FLAG_REV 0x2
#define FLAG_ODD 0x4
#define FLAG_BIN 0x8
#define FLAG_TREE 0x10
#define FLAG_CUTN 0x20
#define FLAG_DROPN 0x40
#define FLAG_RLD 0x80

static void insert1(int flag, int l, uint8_t *s, bprope6_t *bpr, bcr_t *bcr)
{
	int i;
	if ((flag & FLAG_ODD) && (l&1) == 0) { // then check reverse complement
		for (i = 0; i < l>>1; ++i) // is the reverse complement is identical to itself?
			if (s[i] + s[l-1-i] != 5) break;
		if (i == l>>1) --l; // if so, trim 1bp from the end
	}
	if (flag & FLAG_FOR) {
		if (bpr) bpr_insert_string(bpr, l, s);
		if (bcr) bcr_append(bcr, l, s);
	}
	if (flag & FLAG_REV) {
		for (i = 0; i < l>>1; ++i) { // reverse complement
			int tmp = s[l-1-i];
			tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
			s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
			s[i] = tmp;
		}
		if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		if (bpr) bpr_insert_string(bpr, l, s);
		if (bcr) bcr_append(bcr, l, s);
	}
}

int main_ropebwt(int argc, char *argv[])
{
	bprope6_t *bpr = 0;
	bcr_t *bcr = 0;
	gzFile fp;
	FILE *out = stdout;
	char *tmpfn = 0, *outfn = 0;
	kseq_t *ks;
	enum algo_e algo = BPR;
	int c, max_runs = 512, max_nodes = 64, min_len = 1, is_short = 0;
	int bcr_flag = 0, flag = FLAG_FOR | FLAG_REV | FLAG_ODD;

	while ((c = getopt(argc, argv, "TFRODbNo:r:n:ta:f:v:sBl:dS")) >= 0)
		if (c == 'a') {
			if (strcmp(optarg, "bpr") == 0) algo = BPR;
			else if (strcmp(optarg, "bcr") == 0) algo = BCR;
			else fprintf(stderr, "[W::%s] available algorithms: bpr or bcr; default to bpr\n", __func__);
		} else if (c == 'o') out = fopen(optarg, "wb"), outfn = optarg;
		else if (c == 'F') flag &= ~FLAG_FOR;
		else if (c == 'R') flag &= ~FLAG_REV;
		else if (c == 'O') flag &= ~FLAG_ODD;
		else if (c == 'T') flag |= FLAG_TREE;
		else if (c == 'b') flag |= FLAG_BIN;
		else if (c == 'N') flag |= FLAG_CUTN;
		else if (c == 'D') flag |= FLAG_DROPN;
		else if (c == 'd') flag |= FLAG_RLD;
		else if (c == 't') bcr_flag |= BCR_F_THR;
		else if (c == 's') bcr_flag |= BCR_F_RLO;
		else if (c == 'B') bcr_flag |= BCR_F_BPR;
		else if (c == 'l') min_len = atoi(optarg);
		else if (c == 'r') max_runs = atoi(optarg);
		else if (c == 'n') max_nodes= atoi(optarg);
		else if (c == 'f') tmpfn = optarg;
		else if (c == 'v') bcr_verbose = atoi(optarg);
		else if (c == 'S') is_short = 1;
	
	if (is_short) {
		algo = BCR;
		flag = FLAG_FOR | FLAG_REV | FLAG_ODD | FLAG_RLD;// | FLAG_CUTN;
		bcr_flag = BCR_F_THR | BCR_F_RLO;
		min_len = 25;
		bcr_verbose = 3;
	}

	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   ropebwt [options] <in.fq.gz>\n\n");
		fprintf(stderr, "Options: -a STR     algorithm: bpr or bcr [bpr]\n");
		fprintf(stderr, "         -S         build FM-index for short sequences (-Ndstabcr -l25 -v3)\n");
		fprintf(stderr, "         -r INT     max number of runs in leaves (bpr only) [%d]\n", max_runs);
		fprintf(stderr, "         -n INT     max number children per internal node (bpr only) [%d]\n", max_nodes);
		fprintf(stderr, "         -o FILE    output file [stdout]\n");
		fprintf(stderr, "         -f FILE    temporary sequence file name (bcr only) [null]\n");
		fprintf(stderr, "         -v INT     verbose level (bcr only) [%d]\n", bcr_verbose);
		fprintf(stderr, "         -l INT     skip sequences/fragments shorter than INT bp [%d]\n", min_len);
		fprintf(stderr, "         -b         binary output (5+3 runs starting after 4 bytes)\n");
		fprintf(stderr, "         -d         fermi's RLD output\n");
		fprintf(stderr, "         -s         sort reads into RLO order (bcr only)\n");
		fprintf(stderr, "         -B         use rope to keep partial BWT (bcr only)\n");
		fprintf(stderr, "         -t         enable threading (bcr only)\n");
		fprintf(stderr, "         -F         skip forward strand\n");
		fprintf(stderr, "         -R         skip reverse strand\n");
		fprintf(stderr, "         -D         drop reads containing ambiguous bases\n");
		fprintf(stderr, "         -N         cut at ambiguous bases\n");
		fprintf(stderr, "         -O         suppress end trimming when forward==reverse\n");
		fprintf(stderr, "         -T         print the tree stdout (bpr only)\n\n");
		return 1;
	}

	if (algo == BCR) bcr = bcr_init();
	else if (algo == BPR) bpr = bpr_init(max_nodes, max_runs);
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "rb") : gzdopen(fileno(stdin), "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int j;
		uint8_t *t = (uint8_t*)ks->seq.s;
		if (ks->seq.l < min_len) continue;
		for (j = 0; j < ks->seq.l; ++j)
			t[j] = t[j] < 128? seq_nt6_table[t[j]] : 5;
		if (flag & FLAG_CUTN) { // cut at ambiguous bases
			int l;
			uint8_t *s;
			for (j = l = 0, s = t; j < ks->seq.l; ++j) {
				if (t[j] == 5) {
					if (l >= min_len) insert1(flag, l, s, bpr, bcr);
					s = t + j + 1; l = 0;
				} else ++l;
			}
			if (l >= min_len) insert1(flag, l, s, bpr, bcr);
		} else {
			if (flag & FLAG_DROPN) {
				for (j = 0; j < ks->seq.l && t[j] < 5; ++j);
				if (j != ks->seq.l) continue;
			}
			insert1(flag, ks->seq.l, t, bpr, bcr);
		}
	}
	kseq_destroy(ks);
	gzclose(fp);

#define print_bwt(itr_t, itr_set, itr_next_f, is_rld, is_bin, fp) do { \
		itr_t *itr; \
		const uint8_t *s; \
		int i, j, l; \
		itr = (itr_set); \
		if (is_rld) { \
			rld_t *e; \
			rlditr_t ri; \
			e = rld_init(6, 3); \
			rld_itr_init(e, &ri, 0); \
			while ((s = itr_next_f(itr, &l)) != 0) \
				for (i = 0; i < l; ++i) \
					rld_enc(e, &ri, s[i]>>3, s[i]&7); \
			rld_enc_finish(e, &ri); \
			rld_dump(e, outfn? outfn : "-"); \
			rld_destroy(e); \
		} else if (is_bin) { \
			fwrite("RLE\6", 4, 1, fp); \
			while ((s = itr_next_f(itr, &l)) != 0) \
				fwrite(s, 1, l, fp); \
		} else { \
			while ((s = itr_next_f(itr, &l)) != 0) \
				for (i = 0; i < l; ++i) \
					for (j = 0; j < s[i]>>3; ++j) \
						fputc("$ACGTN"[s[i]&7], fp); \
			fputc('\n', fp); \
		} \
		free(itr); \
	} while (0)

	if (bpr) {
		print_bwt(bpriter_t, bpr_iter_init(bpr), bpr_iter_next, flag&FLAG_RLD, flag&FLAG_BIN, out);
		if (flag&FLAG_TREE) bpr_print(bpr);
		bpr_destroy(bpr);
	}
	if (bcr) {
		bcr_build(bcr, bcr_flag, tmpfn);
		print_bwt(bcritr_t, bcr_itr_init(bcr), bcr_itr_next, flag&FLAG_RLD, flag&FLAG_BIN, out);
		bcr_destroy(bcr);
	}
	fclose(out);
	return 0;
}
