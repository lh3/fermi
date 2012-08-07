/* The MIT License

   Copyright (c) 2012 Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef BPROPE6_H
#define BPROPE6_H

#include <stdint.h>

struct bprope6_s;
typedef struct bprope6_s bprope6_t;

struct bpriter_s;
typedef struct bpriter_s bpriter_t;

#ifdef __cplusplus
extern "C" {
#endif

	// return a rope object; $max_nodes: #nodes in an internal bucket; $max_runs: #runs in a leaf
	bprope6_t *bpr_init(int max_nodes, int max_runs);
	// deallocate $rope
	void bpr_destroy(bprope6_t *rope);
	// insert $a after $x symbols in $rope and returns the position of the next insertion
	int64_t bpr_insert_symbol(bprope6_t *rope, int a, int64_t x);
	// insert a string $str of length $l to $rope; NB: a different input order results in a different rope
	void bpr_insert_string(bprope6_t *rope, int l, const uint8_t *str);
	// print the underlying B+ tree of $rope; for debugging only
	void bpr_print(const bprope6_t *rope);
	// ordered iterator
	bpriter_t *bpr_iter_init(const bprope6_t *rope); // to free, simply call free()
	// get the next leaf; $n is set to the number of runs in the leaf
	const uint8_t *bpr_iter_next(bpriter_t *iter, int *n);
	// memory used by the rope
	int64_t bpr_mem(bprope6_t *rope);

#ifdef __cplusplus
}
#endif

#endif
