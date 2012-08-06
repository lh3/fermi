#include <assert.h>
#include <pthread.h>
#include "rld.h"
#include "fermi.h"
#include "utils.h"

static inline void set_bit(uint64_t *bits, uint64_t k)
{
	uint64_t *p = &bits[k>>6];
	k = 1ull << (k&0x3f);
	__sync_or_and_fetch(p, k);
}

static void set_bits(const rld_t *e, const uint64_t *sub, uint64_t *bits, int start, int step)
{
	uint64_t i, k, *ok;
	int c;
	ok = alloca(8 * e->asize);
	for (i = start; i < e->mcnt[1]; i += step) {
		if ((sub[i>>6]>>(i&0x3f)&1) == 0) continue;
		for (k = i;;) {
			set_bit(bits, k);
			c = rld_rank1a(e, k, ok);
			if (c == 0) break;
			else k = e->cnt[c] + ok[c] - 1;
		}
	}
}

static rld_t *gen_idx(rld_t *e0, uint64_t *bits, int is_comp)
{
	int c = 0, c0 = -1;
	int64_t l, i, k = 0, len = 0;
	rld_t *e;
	rlditr_t ritr, witr;

	e = rld_init(e0->asize, e0->sbits);
	rld_itr_init(e, &witr, 0);
	rld_itr_init(e0, &ritr, 0);
	while ((l = rld_dec(e0, &ritr, &c, 1)) >= 0) {
		for (i = 0; i < l; ++i, ++k) {
			if ((bits[k>>6]>>(k&0x3f)&1) == !is_comp) {
				if (c != c0) {
					if (len) rld_enc(e, &witr, len, c0);
					c0 = c, len = 1;
				} else ++len;
			}
		}
	}
	if (len) rld_enc(e, &witr, len, c0);
	assert(k == e0->mcnt[0]);
	rld_destroy(e0);
	rld_enc_finish(e, &witr);
	return e;
}

typedef struct {
	int start, step;
	const rld_t *e;
	const uint64_t *sub;
	uint64_t *bits;
} worker_t;

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	set_bits(w->e, w->sub, w->bits, w->start, w->step);
	return 0;
}

rld_t *fm_sub(rld_t *e, const uint64_t *sub, int n_threads, int is_comp)
{
	uint64_t *bits;
	worker_t *w;
	pthread_t *tid;
	pthread_attr_t attr;
	int i, j;
	rld_t *r;

	if (n_threads < 1) n_threads = 1;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	w = (worker_t*)calloc(n_threads, sizeof(worker_t));
	bits = xcalloc((e->mcnt[0] + 63) / 64, 8);
	for (i = 0; i < n_threads; ++i) {
		w[i].e = e, w[i].sub = sub, w[i].bits = bits;
		w[i].start = i, w[i].step = n_threads;
	}
	for (j = 0; j < n_threads; ++j) pthread_create(&tid[j], &attr, worker, w + j);
	for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
	free(tid); free(w);

	r = gen_idx(e, bits, is_comp);
	free(bits);
	return r;
}
