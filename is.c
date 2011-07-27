/*
 * Copyright (c) 2008 Yuta Mori All Rights Reserved.
 *               2011 Heng Li
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdlib.h>
#include <limits.h>

typedef unsigned char ubyte_t;
/* T is of type "const unsigned char*". If T[i] is a sentinel, chr(i) takes a negative value */
#define chr(i) (cs == sizeof(int) ? ((const int *)T)[i] : (T[i]? (int)T[i] : i-INT_MAX))

/** Count the occurrences of each symbol */
static void getCounts(const unsigned char *T, int *C, int n, int k, int cs)
{
	int i;
	for (i = 0; i < k; ++i) C[i] = 0;
	for (i = 0; i < n; ++i) {
		int c = chr(i);
		++C[c > 0? c : 0];
	}
}

/**
 * Fine the end of each bucket
 *
 * @param C   occurrences computed by getCounts(); input
 * @param B   start/end of each bucket; output
 * @param k   size of alphabet
 * @param end compute the end of bucket if true; otherwise compute the end
 */
static inline void getBuckets(const int *C, int *B, int k, int end)
{
	int i, sum = 0;
	if (end) for (i = 0; i < k; ++i) sum += C[i], B[i] = sum;
	else for (i = 0; i < k; ++i) sum += C[i], B[i] = sum - C[i];
}

/* compute SA */
static void induceSA(const unsigned char *T, int *SA, int *C, int *B, int n, int k, int cs)
{
	int *b, i, j;
	int  c0, c1;
	/* left-to-right induced sort (for L-type) */
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 0);	/* find starts of buckets */
	c1 = chr(SA[0]); b = SA + B[c1 > 0? c1 : 0];
	for (i = 0; i < n; ++i) {
		j = SA[i], SA[i] = ~j;
		if (0 < j) { // >0 if j-1 is L-type; <0 if S-type; ==0 undefined
			--j;
			if ((c0 = chr(j)) != c1) {
				B[c1 > 0? c1 : 0] = b - SA;
				c1 = c0;
				b = SA + B[c1 > 0? c1 : 0];
			}
			*b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
		}
	}
	/* right-to-left induced sort (for S-type) */
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	/* find ends of buckets */
	--B[0]; /* This line to deal with the last sentinel. */
	for (i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
		if (0 < (j = SA[i])) { // the prefix is S-type
			--j;
			if ((c0 = chr(j)) != c1) {
				B[c1 > 0? c1 : 0] = b - SA;
				c1 = c0;
				b = SA + B[c1 > 0? c1 : 0];
			}
			*--b = ((j == 0) || (chr(j - 1) > c1)) ? ~j : j;
		} else SA[i] = ~j; // if L-type, change the sign
	}
}

/**
 * Recursively compute the suffix array.
 *
 * @param T   input string
 * @param SA  output suffix array
 * @param fs  working space available in SA
 * @param n   length of T
 * @param k   size of the alphabet, which may vary with recursion
 * @param cs  # bytes per element in T which can be 1 (uint8_t) or 4 (int32_t)
 */
static int sais_main(const unsigned char *T, int *SA, int fs, int n, int k, int cs)
{
	int *C, *B;
	int  i, j, c, m, q, qlen, name;
	int  c0, c1;

	/* STAGE I: reduce the problem by at least 1/2 sort all the S-substrings */
	if (k <= fs) {
		C = SA + n;
		B = (k <= fs - k) ? C + k : C;
	} else if ((C = (int*)malloc(k * 2 * sizeof(int))) == NULL) return -2;
	else B = C + k;
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	/* find ends of buckets */
	for (i = 0; i < n; ++i) SA[i] = 0;
	/* mark L and S (the t array in Nong et al.), and keep the positions of LMS in the buckets */
	for (i = n - 2, c = 1, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr(i)) < c1 + c) c = 1; /* c1 = chr(i+1); c==1 if in an S run */
		else if (c != 0) SA[--B[c1 > 0? c1 : 0]] = i + 1, c = 0;
	}
	induceSA(T, SA, C, B, n, k, cs);
	if (fs < k) free(C);
	/* pack all the sorted LMS into the first m items of SA 
	   2*m must be not larger than n (see Nong et al. for the proof) */
	for (i = 0, m = 0; i < n; ++i) {
		int p = SA[i];
		if (p == n - 1) SA[m++] = p;
		else if (0 < p && chr(p - 1) > (c0 = chr(p))) {
			for (j = p + 1; j < n && c0 == (c1 = chr(j)); ++j);
			if (j < n && c0 < c1) SA[m++] = p;
		}
	}
	for (i = m; i < n; ++i) SA[i] = 0;	/* init the name array buffer */
	/* store the length of all substrings */
	for (i = n - 2, j = n, c = 1, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr(i)) < c1 + c) c = 1; /* c1 = chr(i+1) */
		else if (c != 0) SA[m + ((i + 1) >> 1)] = j - i - 1, j = i + 1, c = 0;
	}
	/* find the lexicographic names of all substrings */
	for (i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
		int p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
		if (plen == qlen) {
			for (j = 0; (j < plen) && (chr(p + j) == chr(q + j)); j++);
			if (j == plen) diff = 0;
		}
		if (diff != 0) ++name, q = p, qlen = plen;
		SA[m + (p >> 1)] = name;
	}

	/* STAGE II: solve the reduced problem; recurse if names are not yet unique */
	if (name < m) {
		int *RA = SA + n + fs - m;
		for (i = n - 1, j = m - 1; m <= i; --i) {
			if (SA[i] != 0) RA[j--] = SA[i] - 1;
		}
		if (sais_main((unsigned char *)RA, SA, fs + n - m * 2, m, name, sizeof(int)) != 0) return -2;
		for (i = n - 2, j = m - 1, c = 1, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
			if ((c0 = chr(i)) < c1 + c) c = 1;
			else if (c != 0) RA[j--] = i + 1, c = 0; /* get p1 */
		}
		for (i = 0; i < m; ++i) SA[i] = RA[SA[i]]; /* get index */
	}

	/* STAGE III: induce the result for the original problem */
	if (k <= fs) {
		C = SA + n;
		B = (k <= fs - k) ? C + k : C;
	} else if ((C = (int *)malloc(k * 2 * sizeof(int))) == NULL) return -2;
	else B = C + k;
	/* put all LMS characters into their buckets */
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	/* find ends of buckets */
	for (i = m; i < n; ++i) SA[i] = 0; /* init SA[m..n-1] */
	for (i = m - 1; 0 <= i; --i) {
		j = SA[i], SA[i] = 0;
		c = chr(j);
		SA[--B[c > 0? c : 0]] = j;
	}
	induceSA(T, SA, C, B, n, k, cs);
	if (fs < k) free(C);
	return 0;
}

/**
 * Constructs the suffix array of a given string.
 *
 * @param T[0..n-1]  NULL terminated input string
 * @param SA[0..n-1] Output suffix array
 * @param n          The length of the given string.
 * @return 0         If no error occurred
 */
int is_sa(const ubyte_t *T, int *SA, int n)
{
	if ((T == NULL) || (SA == NULL) || (n <= 0)) return -1;
	return sais_main(T, SA, 0, n, 256, 1);
}
