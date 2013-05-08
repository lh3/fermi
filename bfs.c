#include <string.h>
#include "fermi.h"
#include "rld.h"

typedef struct {
	uint64_t k:63, f:1;
	uint64_t l:56, d:8;
} fm6_bfs1_t;

typedef struct {
	int bits;
	int64_t front, count;
	fm6_bfs1_t *a;
} queue_t;

static queue_t *qu_init()
{
	queue_t *q;
	q = calloc(1, sizeof(queue_t));
	q->bits = 2;
	q->a = malloc((1<<q->bits) * sizeof(fm6_bfs1_t));
	return q;
}

static inline fm6_bfs1_t *qu_enq(queue_t *q) // en-queue
{
	if (q->count == 1LL<<q->bits) { // re-queue
		++q->bits;
		q->a = realloc(q->a, (1LL<<q->bits) * sizeof(fm6_bfs1_t));
		memmove(q->a + q->count, q->a, q->front * sizeof(fm6_bfs1_t));
		memmove(q->a, q->a + q->front, q->count * sizeof(fm6_bfs1_t));
		q->front = 0;
	}
	return &q->a[((q->count++) + q->front) & ((1LL<<q->bits) - 1)];
}

static inline void qu_deq(queue_t *q, fm6_bfs1_t *d) // de-queue
{
	if (q->count == 0) return;
	*d = q->a[q->front++];
	q->front &= (1LL<<q->bits) - 1;
	--q->count;
}

void fm6_bfs_traverse(const rld_t *e, uint64_t k0, uint64_t l0)
{
	uint64_t *bl, *br;
	queue_t *queue;
	fm6_bfs1_t *p;

	queue = qu_init();
	bl = calloc((e->mcnt[0] + 64) / 64, 8);
	br = calloc((e->mcnt[0] + 64) / 64, 8);
	p = qu_enq(queue);
	p->k = k0; p->l = l0; p->d = 0;
	while (queue->count) {
		uint64_t yl, yr, *pl, *pr, zl, zr;
		uint64_t ok[6], ol[6];
		fm6_bfs1_t x;
		int c;

		qu_deq(queue, &x);
//		printf("%ld\t%ld\t%ld\n", (long)x.d, (long)x.k, (long)x.l);
		rld_rank2a(e, (int64_t)x.k - 1, x.l - 1, ok, ol);
		for (c = 1; c < e->asize; ++c) {
			if (ok[c] == ol[c]) continue;
			x.k = e->cnt[c] + ok[c];
			x.l = e->cnt[c] + ol[c];
			pl = bl + (x.k>>6); zl = 1ULL << (x.k&0x3f);
			pr = br + (x.l>>6); zr = 1ULL << (x.l&0x3f);
			yl = __sync_fetch_and_or(pl, zl);
			yr = __sync_fetch_and_or(pr, zr);
			if ((yl&zl) && (yr&zr)) continue;
			p = qu_enq(queue);
			p->k = e->cnt[c] + ok[c];
			p->l = e->cnt[c] + ol[c];
			p->d = x.d < 255? x.d + 1 : 255;
		}
	}
	free(bl); free(br);
	free(queue->a); free(queue);
}
