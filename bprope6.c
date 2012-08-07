#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "bprope6.h"

#define MP_CHUNK_SIZE 0x100000 // 1MB per chunk

typedef struct { // memory pool for fast and compact memory allocation (no free)
	int size, i, n_elems;
	int64_t top, max;
	uint8_t **mem;
} mempool_t;

static mempool_t *mp_init(int size)
{
	mempool_t *mp;
	mp = calloc(1, sizeof(mempool_t));
	mp->size = size;
	mp->i = mp->n_elems = MP_CHUNK_SIZE / size;
	mp->top = -1;
	return mp;
}

static void mp_destroy(mempool_t *mp)
{
	int64_t i;
	for (i = 0; i <= mp->top; ++i) free(mp->mem[i]);
	free(mp->mem); free(mp);
}

static inline void *mp_alloc(mempool_t *mp)
{
	if (mp->i == mp->n_elems) {
		if (++mp->top == mp->max) {
			mp->max = mp->max? mp->max<<1 : 1;
			mp->mem = realloc(mp->mem, sizeof(void*) * mp->max);
		}
		mp->mem[mp->top] = calloc(mp->n_elems, mp->size);
		mp->i = 0;
	}
	return mp->mem[mp->top] + (mp->i++) * mp->size;
}

static int insert_to_leaf(uint8_t *p, int a, int x, int len, uint64_t c[6])
{ // insert $a after $x symbols in $p; IMPORTANT: the first 4 bytes of $p gives the length of the string
#define MAX_RUNLEN 31
#define _insert_after(_n, _s, _i, _b) if ((_i) + 1 != (_n)) memmove(_s+(_i)+2, _s+(_i)+1, (_n)-(_i)-1); _s[(_i)+1] = (_b); ++(_n)

	int r[6], i, l = 0, n = *(int32_t*)p;
	uint8_t *s = p + 4;
	if (n == 0) { // if $s is empty, that is easy
		s[n++] = 1<<3 | a;
		*(int32_t*)p = n;
		return 0;
	}
	if (x < len>>1) { // forwardly search for the run to insert
		for (i = 0; i < 6; ++i) r[i] = 0;
		do {
			l += *s>>3;
			r[*s&7] += *s>>3;
			++s;
		} while (l < x);
	} else { // backwardly search for the run to insert; this block has exactly the same functionality as the above
		for (i = 0; i < 6; ++i) r[i] = c[i];
		l = len, s += n;
		do {
			--s;
			l -= *s>>3;
			r[*s&7] -= *s>>3;
		} while (l >= x);
		l += *s>>3; r[*s&7] += *s>>3; ++s;
	}
	i = s - p - 4; s = p + 4;
	assert(i <= n);
	r[s[--i]&7] -= l - x; // $i now points to the left-most run where $a can be inserted
	if (l == x && i != n - 1 && (s[i+1]&7) == a) ++i; // if insert to the end of $i, check if we'd better to the start of ($i+1)
	if ((s[i]&7) == a) { // insert to a long $a run
		if (s[i]>>3 == MAX_RUNLEN) { // the run is full
			for (++i; i != n && (s[i]&7) == a; ++i); // find the end of the long run
			--i;
			if (s[i]>>3 == MAX_RUNLEN) { // then we have to add one run
				_insert_after(n, s, i, 1<<3|a);
			} else s[i] += 1<<3;
		} else s[i] += 1<<3;
	} else if (l == x) { // insert to the end of run; in this case, neither this and the next run is $a
		_insert_after(n, s, i, 1<<3 | a);
	} else if (i != n - 1 && (s[i]&7) == (s[i+1]&7)) { // insert to a long non-$a run
		int rest = l - x, c = s[i]&7;
		s[i] -= rest<<3;
		_insert_after(n, s, i, 1<<3 | a);
		for (i += 2; i != n && (s[i]&7) == c; ++i); // find the end of the long run
		--i;
		if ((s[i]>>3) + rest > MAX_RUNLEN) { // we cannot put $rest to $s[$i]
			rest = (s[i]>>3) + rest - MAX_RUNLEN;
			s[i] = MAX_RUNLEN<<3 | (s[i]&7);
			_insert_after(n, s, i, rest<<3 | c);
		} else s[i] += rest<<3;
	} else { // insert to a short run
		memmove(s + i + 3, s + i + 1, n - i - 1);
		s[i]  -= (l-x)<<3;
		s[i+1] = 1<<3 | a;
		s[i+2] = (l-x)<<3 | (s[i]&7);
		n += 2;
	}
	*(int32_t*)p = n;
	return r[a];
}

typedef struct bpr_node_s {
	struct bpr_node_s *p; // child; at the bottom level, $p points to a string with the first 4 bytes giving the number of runs (#runs)
	uint64_t l:54, n:9, is_bottom:1; // $n and $is_bottom are only set for the first node in a bucket
	uint64_t c[6]; // marginal counts
} node_t;

struct bprope6_s {
	int max_nodes, max_runs; // both MUST BE even numbers
	uint64_t c[6]; // marginal counts
	node_t *root;
	mempool_t *node, *leaf;
};

static void print_node(const node_t *p) // recursively print the B+ rope in the Newick format
{
	if (p->is_bottom) {
		int i, j, k;
		putchar('(');
		for (i = 0; i < p->n; ++i) {
			uint8_t *q = (uint8_t*)p[i].p;
			int n = *(int32_t*)q;
			if (i) putchar(',');
			for (j = 0, q += 4; j < n; ++j)
				for (k = 0; k < q[j]>>3; ++k)
					putchar("$ACGTN"[q[j]&7]);
		}
		putchar(')');
	} else {
		int i;
		putchar('(');
		for (i = 0; i < p->n; ++i) {
			if (i) putchar(',');
			print_node(p[i].p);
		}
		putchar(')');
	}
}

void bpr_print(const bprope6_t *rope) { print_node(rope->root); putchar('\n'); }

static inline node_t *split_node(bprope6_t *rope, node_t *u, node_t *v)
{ // split $v's child. $u is the first node in the bucket. $v and $u are in the same bucket. IMPORTANT: there is always enough room in $u
	int j, i = v - u;
	node_t *w; // $w is the sibling of $v
	if (u == 0) { // only happens at the root; add a new root
		u = v = mp_alloc(rope->node);
		v->n = 1; v->p = rope->root; // the new root has the old root as the only child
		memcpy(v->c, rope->c, 48);
		for (j = 0; j < 6; ++j) v->l += v->c[j];
		rope->root = v;
	}
	if (i != u->n - 1) // then make room for a new node
		memmove(v + 2, v + 1, sizeof(node_t) * (u->n - i - 1));
	++u->n; w = v + 1;
	memset(w, 0, sizeof(node_t));
	w->p = mp_alloc(u->is_bottom? rope->leaf : rope->node);
	if (u->is_bottom) { // we are at the bottom level; $v->p is a string instead of a node
		uint8_t *p = (uint8_t*)v->p, *q = (uint8_t*)w->p;
		int32_t *np = (int32_t*)p, *nq = (int32_t*)q; // the first 4 bytes give the number of runs (#runs)
		*nq = *np - (rope->max_runs>>1); *np -= *nq; // compute the new #runs
		memcpy(q + 4, p + 4 + *np, *nq);
		for (i = 0, q += 4; i < *nq; ++i) // compute the count arrays in $w, from $q
			w->c[q[i]&7] += q[i]>>3;
	} else { // $v->p is a node, not a string
		node_t *p = v->p, *q = w->p; // $v and $w are siblings and thus $p and $q are cousins
		p->n -= rope->max_nodes>>1;
		memcpy(q, p + p->n, sizeof(node_t) * (rope->max_nodes>>1));
		q->n = rope->max_nodes>>1; // NB: this line must below memcpy() as $q->n and $q->is_bottom are modified by memcpy()
		q->is_bottom = p->is_bottom;
		for (i = 0; i < q->n; ++i)
			for (j = 0; j < 6; ++j)
				w->c[j] += q[i].c[j];
	}
	for (j = 0; j < 6; ++j) // compute $w->l and update $v->c
		w->l += w->c[j], v->c[j] -= w->c[j];
	v->l -= w->l; // update $v->c
	return v;
}

int64_t bpr_insert_symbol(bprope6_t *rope, int a, int64_t x)
{ // insert $a after $x symbols in $rope and the returns the position of the next insertion
	node_t *u = 0, *v = 0, *p = rope->root; // $v is the parent of $p; $u and $v are at the same level and $u is the first node in the bucket
	int64_t y = 0, z;
	int i;
	for (i = 0, z = 0; i < a; ++i) z += rope->c[i];
	do { // top-down update. Searching and node splitting are done together in one pass.
		if (p->n == rope->max_nodes) { // node is full; split
			v = split_node(rope, u, v); // $v points to the parent of $p; when a new root is added, $v points to the root
			if (y + v->l < x) // if $v is not long enough after the split, we need to move both $p and its parent $v
				y += v->l, z += v->c[a], ++v, p = v->p;
		}
		u = p;
		if (v && x - y > v->l>>1) { // then search backwardly for the right node to descend
			p += p->n - 1; y += v->l; z += v->c[a];
			for (; y >= x; --p) y -= p->l, z -= p->c[a];
			++p;
		} else for (; y + p->l < x; ++p) y += p->l, z += p->c[a]; // then search forwardly
		assert(p - u < u->n);
		if (v) ++v->c[a], ++v->l; // we should not change p->c[a] because this may cause troubles when p's child is split
		v = p; p = p->p; // descend
	} while (!u->is_bottom);
	++rope->c[a]; // $rope->c should be updated after the loop as adding a new root needs the old $rope->c counts
	z += insert_to_leaf((uint8_t*)p, a, x - y, v->l, v->c) + 1;
	++v->c[a]; ++v->l; // this should be below insert_to_leaf(); otherwise insert_to_leaf() will not work
	if (*(uint32_t*)p + 2 > rope->max_runs) split_node(rope, u, v);
	return z;
}

void bpr_insert_string(bprope6_t *rope, int l, const uint8_t *str)
{
	uint64_t x = rope->c[0];
	for (--l; l >= 0; --l)
		x = bpr_insert_symbol(rope, str[l], x);
	bpr_insert_symbol(rope, 0, x);
}

bprope6_t *bpr_init(int max_nodes, int max_runs)
{
	bprope6_t *rope;
	rope = calloc(1, sizeof(bprope6_t));
	if (max_runs < 8) max_runs = 8;
	rope->max_nodes= (max_nodes+ 1)>>1<<1;
	rope->max_runs = ((max_runs + 1)>>1<<1) - 4; // -4 to make room for the 4-byte integer keeping #runs
	rope->node = mp_init(sizeof(node_t) * rope->max_nodes);
	rope->leaf = mp_init(rope->max_runs + 4); // +4 to include the number of runs
	rope->root = mp_alloc(rope->node);
	rope->root->n = 1;
	rope->root->is_bottom = 1;
	rope->root->p = mp_alloc(rope->leaf);
	return rope;
}

void bpr_destroy(bprope6_t *rope)
{
	mp_destroy(rope->node);
	mp_destroy(rope->leaf);
	free(rope);
}

int64_t bpr_mem(bprope6_t *rope)
{
	return (rope->leaf->top + rope->node->top + 2) * MP_CHUNK_SIZE + (rope->leaf->max + rope->node->max) * sizeof(void*)
		+ sizeof(mempool_t) * 2 + sizeof(bprope6_t) + 10 * sizeof(void*);
}

struct bpriter_s {
	const bprope6_t *rope;
	const node_t *pa[80];
	int k, ia[80];
};

bpriter_t *bpr_iter_init(const bprope6_t *rope)
{
	bpriter_t *i;
	i = calloc(1, sizeof(bpriter_t));
	i->rope = rope;
	for (i->pa[i->k] = rope->root; !i->pa[i->k]->is_bottom;) // descend to the leftmost leaf
		++i->k, i->pa[i->k] = i->pa[i->k - 1]->p;
	return i;
}

const uint8_t *bpr_iter_next(bpriter_t *i, int *n)
{
	const uint8_t *ret;
	assert(i->k < 80); // a B+ tree should not be that tall
	if (i->k < 0) return 0;
	*n = *(int32_t*)i->pa[i->k][i->ia[i->k]].p;
	ret = (uint8_t*)i->pa[i->k][i->ia[i->k]].p + 4;
	while (i->k >= 0 && ++i->ia[i->k] == i->pa[i->k]->n) i->ia[i->k--] = 0; // backtracking
	if (i->k >= 0)
		while (!i->pa[i->k]->is_bottom) // descend to the leftmost leaf
			++i->k, i->pa[i->k] = i->pa[i->k - 1][i->ia[i->k - 1]].p;
	return ret;
}
