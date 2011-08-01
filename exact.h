#ifndef BWA_SEARCH_H
#define BWA_SEARCH_H

#include "rld.h"
#include "kstring.h"

#ifdef __cplusplus
extern "C" {
#endif

	uint64_t fm_backward_search(rld_t *e, const rldidx_t *r, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end);
	void fm_retrieve(rld_t *e, const rldidx_t *r, uint64_t x, kstring_t *s);

#ifdef __cplusplus
}
#endif

#endif
