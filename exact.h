#ifndef BWA_SEARCH_H
#define BWA_SEARCH_H

#include "rld.h"
#include "kstring.h"

#ifdef __cplusplus
extern "C" {
#endif

	uint64_t fm_backward_search(const rld_t *e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end);
	void fm_retrieve(const rld_t *e, uint64_t x, kstring_t *s);
	void fm6_retrieve(const rld_t *e, uint64_t x, kstring_t *s);

#ifdef __cplusplus
}
#endif

#endif
