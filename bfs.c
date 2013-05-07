#include <string.h>
#include "fermi.h"
#include "rld.h"

typedef struct {
	int bits;
	int64_t front, count;
	ku128_t *a;
} queue_t;

static queue_t *qu_init()
{
	queue_t *q;
	q = calloc(1, sizeof(queue_t));
	q->bits = 2;
	q->a = malloc((1<<q->bits) * sizeof(ku128_t));
	return q;
}

static inline ku128_t *qu_enq(queue_t *q) // en-queue
{
	if (q->count == 1LL<<q->bits) { // re-queue
		++q->bits;
		q->a = realloc(q->a, (1LL<<q->bits) * sizeof(ku128_t));
		memmove(q->a + q->count, q->a, q->front * sizeof(ku128_t));
		memmove(q->a, q->a + q->front, q->count * sizeof(ku128_t));
		q->front = 0;
	}
	return &q->a[((q->count++) + q->front) & ((1LL<<q->bits) - 1)];
}

static inline void qu_deq(queue_t *q, ku128_t *d) // de-queue
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
	ku128_t *p;

	queue = qu_init();
	list = kl_init(intv);
	bl = calloc((e->mcnt[0] + 64) / 64, 8);
	br = calloc((e->mcnt[0] + 64) / 64, 8);
	p = qu_enq(queue);
	p->x = k0; p->y = l0;
	while (queue->count) {
		uint64_t yl, yr, *pl, *pr, zl, zr;
		uint64_t ok[6], ol[6];
		ku128_t x;
		int c;

		qu_deq(queue, &x);
		printf("[%ld,%ld)\n", (long)x.x, (long)x.y);
		rld_rank2a(e, x.x - 1, x.y - 1, ok, ol);
		for (c = 1; c < e->asize; ++c) {
			if (ok[c] == ol[c]) continue;
			x.x = e->cnt[c] + ok[c];
			x.y = e->cnt[c] + ol[c];
			pl = bl + (x.x>>6); zl = 1ULL << (x.x&0x3f);
			pr = br + (x.y>>6); zr = 1ULL << (x.y&0x3f);
			yl = __sync_fetch_and_or(pl, zl);
			yr = __sync_fetch_and_or(pr, zr);
			if ((yl&zl) && (yr&zr)) continue;
			p = qu_enq(queue);
			p->x = e->cnt[c] + ok[c];
			p->y = e->cnt[c] + ol[c];
		}
	}
	free(bl); free(br);
	free(queue->a); free(queue);
}
