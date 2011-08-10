#ifndef FERMI_H
#define FERMI_H

#include <stdint.h>

#define FERMI_VERSION "0.0-dev (r102)"

typedef struct {
	uint64_t x[3];
} fmintv_t;

struct __rld_t;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	uint32_t l, m;
	char *s;
} kstring_t;
#endif

/* complement of a nucleotide */
#define fm6_comp(a) ((a) >= 1 && (a) <= 4? 5 - (a) : (a))

#ifdef __cplusplus
extern "C" {
#endif

	uint64_t fm_backward_search(const struct __rld_t *e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end);
	void fm_retrieve(const struct __rld_t *e, uint64_t x, kstring_t *s);
	void fm6_retrieve(const struct __rld_t *e, uint64_t x, kstring_t *s);

	int fm6_extend(const struct __rld_t *e, const fmintv_t *ik, fmintv_t ok[6], int is_back);

	int fm6_search_forward_overlap(const struct __rld_t *e, int min, int len, const uint8_t *seq);

	struct __rld_t *fm_merge0(const struct __rld_t *e0, const struct __rld_t *e1);

#ifdef __cplusplus
}
#endif

#endif
