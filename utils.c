#include <sys/resource.h>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

// 1: error; 2: warning; 3: message; 4: progress; 5: debugging; >=10: pure debugging
int fm_verbose = 4;

void *xmalloc_core(size_t s, const char *func)
{
	unsigned char *x;
	x = (unsigned char*)malloc(s);
	if (x == 0) {
		fprintf(stderr, "[E::%s] Fail to allocate %ld bytes of memory.\n", func, s);
		return 0;
	}
	return (void*)x;
}

void *xcalloc_core(size_t n, size_t s, const char *func)
{
	unsigned char *x;
	x = (unsigned char*)calloc(n, s);
	if (x == 0) {
		fprintf(stderr, "[E::%s] Fail to allocate %ld bytes of memory.\n", func, s * n);
		return 0;
	}
	return (void*)x;
}

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

/* Comment on the memory usage. We may think getrusage() is the easiest way to
 * get the memory usage of the current process.  However, on Linux, ru_maxrss
 * is not supported by the kernel. We can get the current usage from file
 * "/proc/self/stat" if the "/proc" file system is built in the kernel, but
 * this is a snapshot instead of the maximum memory usage. We have to monitor
 * the "stat" file constantly to get the maximum resident memory (RSS). On
 * Linux, perhaps the easiest way to get the maximum RSS is to use my "runit".
 * Even "/usr/bin/time" may not always work on all versions/distributions.
 *
 * On Mac, getrusage() writes ru_maxrss, but it seems to me that this field
 * tends to be much larger than my expectation. This may be caused by the
 * memory manager (i.e. malloc) instead of the kernel. When I have time,
 * perhaps it is worth investigating an alternative memory manager or
 * controling the major memory allocation by myself, which though could make
 * the code look nasty.
 */
#if defined(__linux__)
double rssmem()
{
	FILE *fp;
	int c, n_spc = 0;
	long mem, page_size;
	page_size = sysconf(_SC_PAGESIZE);
	fp = fopen("/proc/self/stat", "r");
	while ((c = fgetc(fp)) != EOF) {
		if (c == ' ') ++n_spc;
		if (n_spc == 23) break;
	}
	fscanf(fp, "%ld", &mem);
	fclose(fp);
	return mem * page_size / 1024.0 / 1024.0;
}
#elif defined(__APPLE__)
double rssmem()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_maxrss / 1024.0 / 1024.0;
}
#else
double rssmem() { return 0.; }
#endif

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
