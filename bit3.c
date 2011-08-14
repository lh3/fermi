#include <string.h>
#include "bit3.h"

bit3_t *b3_init()
{
	bit3_t *b3;
	int i, j;
	b3 = calloc(1, sizeof(bit3_t));
	b3->n = 1;
	b3->z = malloc(sizeof(void*));
	b3->z[0] = calloc(B3_LSIZE, 2);
	b3->mcnt = calloc(7, 8);
	b3->cnt = calloc(7, 8);
	b3->cnttab = calloc(0x10000, 8);
	for (i = 0; i < 0x10000; ++i) {
		uint64_t x = 0;
		for (j = 0; j <= 12; j += 3)
			x += 1ull<<((i>>j&7)<<3);
		b3->cnttab[i] = x;
	}
	return b3;
}

void b3_destroy(bit3_t *b3)
{
	int i;
	for (i = 0; i < b3->n; ++i) free(b3->z[i]);
	free(b3->z); free(b3->cnt); free(b3->mcnt); free(b3->cnttab); free(b3->frame); free(b3);
}

int64_t b3_enc_finish(bit3_t *b3, b3itr_t *itr)
{
	int i, j, stop = 0;
	uint64_t *q;
	b3_enc(b3, itr, 1, 7);
	while (itr->r || (itr->p - *itr->i) % B3_SFSIZE) b3_enc(b3, itr, 1, 7);
	b3->n_bytes = ((uint64_t)(b3->n - 1) * B3_LSIZE + (itr->p - *itr->i)) * 2;
	b3->n_frames = b3->n_bytes / B3_SSIZE / B3_FSIZE * 6 + b3->n_bytes / B3_SSIZE;
	b3->frame = q = malloc(b3->n_frames * 8);
	itr->i = b3->z; itr->p = *itr->i;
	for (i = 0; i < 7; ++i) *q++ = 0;
	do {
		for (i = 0; i < B3_FSIZE; ++i) {
			uint64_t c = 0;
			for (j = 0; j < B3_SSIZE; ++j)
				c += b3->cnttab[*itr->p++];
			*q++ = c;
			if (c>>56) stop = 1;
			for (j = 0; j < 6; ++j, c >>= 8)
				b3->cnt[j+1] += c&0xff;
			b3->cnt[0] += 3 * B3_SSIZE;
		}
		for (j = 0; j < 6; ++j)
			*q++ = b3->cnt[j+1];
		if (itr->p - *itr->i == B3_LSIZE) itr->p = *++itr->i;
	} while (!stop);
	for (j = 0; j < 7; ++j) b3->mcnt[j] = b3->cnt[j];
	b3->cnt[0] = 0;
	for (j = 2; j < 7; ++j) b3->cnt[j] += b3->cnt[j-1];
	return b3->n_bytes;
}

int b3_dump(const bit3_t *b3, const char *fn)
{
	int64_t k;
	int i;
	FILE *fp;
	fp = strcmp(fn, "-")? fopen(fn, "wb") : stdout;
	fwrite("BIT\3", 1, 4, fp);
	fwrite(b3->mcnt, 8, 7, fp);
	fwrite(&b3->n_bytes, 8, 1, fp);
	for (i = 0, k = b3->n_bytes/2; i < b3->n - 1; ++i, k -= B3_LSIZE)
		fwrite(b3->z[i], 2, B3_LSIZE, fp);
	fwrite(b3->z[i], 2, k, fp);
	fwrite(&b3->n_frames, 8, 1, fp);
	fwrite(&b3->frame, 8, b3->n_frames, fp);
	fclose(fp);
	return 0;
}

bit3_t *b3_restore(FILE *fp)
{
	bit3_t *b3;
	int i;
	int64_t k;
	char magic[5];
	b3 = b3_init();
	fread(magic, 1, 4, fp);
	fread(b3->mcnt, 8, 7, fp);
	for (i = 1; i < 7; ++i) b3->cnt[i] = b3->mcnt[i];
	for (i = 1; i < 7; ++i) b3->cnt[i] += b3->cnt[i-1];
	fread(&b3->n_bytes, 8, 1, fp);
	if (b3->n_bytes / 2 > B3_LSIZE) {
		b3->n = (b3->n_bytes / 2 + B3_LSIZE - 1) / B3_LSIZE;
		b3->z = realloc(b3->z, b3->n * sizeof(void*));
		for (i = 1; i < b3->n; ++i)
			b3->z[i] = calloc(B3_LSIZE, 2);
	}
	for (i = 0, k = b3->n_bytes / 2; i < b3->n - 1; ++i, k -= B3_LSIZE)
		fread(b3->z[i], 2, B3_LSIZE, fp);
	fread(b3->z[i], 2, k, fp);
	fread(&b3->n_frames, 8, 1, fp);
	b3->frame = malloc(b3->n_frames * 8);
	fread(b3->frame, 8, b3->n_frames, fp);
	return b3;
}
