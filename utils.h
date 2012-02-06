#ifndef FERMI_UTILS_H

#include <stdint.h>

#define xmalloc(s) xmalloc_core((s), __func__)
#define xcalloc(n, s) xcalloc_core((n), (s), __func__)

#ifndef KINT_DEF
#define KINT_DEF
typedef struct { uint64_t x, y; } ku128_t;
typedef struct { size_t n, m; uint64_t *a; } ku64_v;
typedef struct { size_t n, m; ku128_t *a; } ku128_v;
#endif

#ifdef __cplusplus
#endif

	double rssmem();
	double cputime();
	double realtime();
	void *xmalloc_core(size_t s, const char *func);
	void *xcalloc_core(size_t n, size_t s, const char *func);

#ifdef __cplusplus
#endif

#endif
