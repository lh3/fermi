#include <sys/resource.h>
#include <sys/time.h>

// 1: error; 2: warning; 3: message
int fm_verbose = 3;

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double rssmem()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_maxrss / 1024.0;
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
