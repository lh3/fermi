#ifndef FERMI_H
#define FERMI_H

#include <stdint.h>

typedef uint64_t fmintv_t[3];

struct __rld_t;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	uint32_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

	uint64_t fm_backward_search(const struct __rld_t *e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end);
	void fm_retrieve(const struct __rld_t *e, uint64_t x, kstring_t *s);
	void fm6_retrieve(const struct __rld_t *e, uint64_t x, kstring_t *s);

	int fm6_extend(const struct __rld_t *e, fmintv_t ik, fmintv_t ok[6], int is_back);

#ifdef __cplusplus
}
#endif

#endif
