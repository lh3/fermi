#include <limits.h>

// modified from suftest.c in libdivsufsort
int sa_check(const unsigned char *T, const int *SA, int n)
{
	#define chr(i) (T[i]? (int)T[i] : i - INT_MAX)
	int C[256];
	int i, p, q, t, c;

	for (i = 0; i < n; ++i)
		if (SA[i] < 0 || n <= SA[i])
			return -2;
	for (i = 1; i < n; ++i)
		if (chr(SA[i - 1]) > chr(SA[i]))
			return -3;
	for (i = 0; i < 256; ++i) C[i] = 0;
	for (i = 0; i < n; ++i) ++C[T[i]];
	for (i = 0, p = 0; i < 256; ++i) t = C[i], C[i] = p, p += t;
	q = C[T[n - 1]];
	++C[T[n - 1]];
	for (i = 0; i < n; ++i) {
		p = SA[i];
		if (0 < p) c = T[--p], t = C[c];
		else c = T[p = n - 1], t = q;
		if (t < 0 || (p != SA[t] && (T[p] || T[SA[t]])))
			return -4;
		if (t != q) {
			++C[c];
			if (n <= C[c] || T[SA[C[c]]] != c)
				C[c] = -1;
		}
	}
	return 0;
	#undef chr
}
