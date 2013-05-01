#ifndef INTHASH_H
#define INTHASH_H

typedef struct {
	volatile int lock;
	uint8_t b_key, b_val, b_tbl, b_idx;
	uint32_t n_elem;
	uint32_t *h;
} inthash_t;

#ifdef __cplusplus
extern "C" {
#endif

	inthash_t *ih_init(int b_key, int b_idx, int b_val);
	uint32_t *ih_put(inthash_t *ih, uint64_t x);
	const uint32_t *ih_get(const inthash_t *ih, uint64_t x);

#ifdef __cplusplus
}
#endif

#define ih_destroy(ih) do { free((ih)->h); free(ih); } while (0)

static inline uint64_t ih_dec(const inthash_t *ih, const uint32_t *p)
{
	uint64_t x;
	x = *p >> ih->b_val & ((1U<<ih->b_idx) - 1);
	x = ((uint64_t)(p - ih->h) - ((x - 1) * x >> 1)) & ((1ULL<<ih->b_tbl) - 1);
	return *p >> (ih->b_idx + ih->b_val) << ih->b_tbl | x;
}

#endif
