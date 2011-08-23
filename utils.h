#ifndef FERMI_UTILS_H

#define xmalloc(s) xmalloc_core((s), __func__)
#define xcalloc(n, s) xcalloc_core((n), (s), __func__)

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
