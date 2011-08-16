#ifndef FERMI_UTILS_H

#define xmalloc(s) xmalloc_core(s, __func__)

#ifdef __cplusplus
#endif

	double cputime();
	double realtime();
	void *xmalloc_core(size_t s, const char *func);

#ifdef __cplusplus
#endif

#endif
