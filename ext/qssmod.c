/* QSufSort.c

   Original source from qsufsort.c

   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

   Modified by Wong Chi-Kwong, 2004

   Changes summary:	- Used long variable and function names
					- Removed global variables
					- Replace pointer references with array references
					- Used insertion sort in place of selection sort and increased insertion sort threshold
					- Reconstructing suffix array from inverse becomes an option
					- Add handling where end-of-text symbol is not necessary < all characters
					- Removed codes for supporting alphabet size > number of characters
  
  No warrenty is given regarding the quality of the modifications.

  Modified by Heng Li, 2011

  Tailor QSufSort for a specific use. The tailored version is NOT able to
  construct the suffix array for ordinary inputs.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define KEY(V, I, p, h)	     ( V[ I[p] + h ] )
#define INSERT_SORT_NUM_ITEM 16

#define min(value1, value2)  ( ((value1) < (value2)) ? (value1) : (value2) )
#define med3(a, b, c)        ( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define swap(a, b, t)        { t = a; a = b; b = t; }

// Static functions
static void QSufSortSortSplit(int* __restrict V, int* __restrict I, const int lowestPos, 
							  const int highestPos, const int numSortedChar);
static int QSufSortChoosePivot(int* __restrict V, int* __restrict I, const int lowestPos, 
							   const int highestPos, const int numSortedChar);
static void QSufSortInsertSortSplit(int* __restrict V, int* __restrict I, const int lowestPos, 
									const int highestPos, const int numSortedChar);

void qsufsort_mod(int *V, int *I, int numChar, int largestInputSymbol, int smallestInputSymbol)
{

	int i, j;
	int s, negatedSortedGroupLength;
	int maxNumInputSymbol;
	int numSortedPos = 1;
   
	maxNumInputSymbol = largestInputSymbol - smallestInputSymbol + 1;

	while ((int)(I[0]) >= -(int)numChar) {
		i = 0;
		negatedSortedGroupLength = 0;
		do {
			s = I[i];
			if (s < 0) {
				i -= s;						/* skip over sorted group.*/
				negatedSortedGroupLength += s;
			} else {
				if (negatedSortedGroupLength) {
					I[i+negatedSortedGroupLength] = negatedSortedGroupLength;	/* combine preceding sorted groups */
					negatedSortedGroupLength = 0;
				}
				j = V[s] + 1;
				QSufSortSortSplit(V, I, i, j - 1, numSortedPos);
				i = j;
			}
		} while (i <= numChar);
		if (negatedSortedGroupLength) {
			/* array ends with a sorted group.*/
			I[i+negatedSortedGroupLength] = negatedSortedGroupLength;	/* combine sorted groups at end of I.*/
		}
		numSortedPos *= 2;	/* double sorted-depth.*/
	}

}

/* Sorting routine called for each unsorted group. Sorts the array of integers
   (suffix numbers) of length n starting at p. The algorithm is a ternary-split
   quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
   Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
   function is based on Program 7.*/
static void QSufSortSortSplit(int* __restrict V, int* __restrict I, const int lowestPos, 
							  const int highestPos, const int numSortedChar) {

	int a, b, c, d;
	int l, m;
	int f, v, s, t;
	int tmp;
	int numItem;

	numItem = highestPos - lowestPos + 1;

	if (numItem <= INSERT_SORT_NUM_ITEM) {
		QSufSortInsertSortSplit(V, I, lowestPos, highestPos, numSortedChar);
		return;
	}

	v = QSufSortChoosePivot(V, I, lowestPos, highestPos, numSortedChar);

	a = b = lowestPos;
	c = d = highestPos;

	while (1) {
		while (c >= b && (f = KEY(V, I, b, numSortedChar)) <= v) {
			if (f == v) {
				swap(I[a], I[b], tmp);
				a++;
			}
			b++;
		}
		while (c >= b && (f = KEY(V, I, c, numSortedChar)) >= v) {
			if (f == v) {
				swap(I[c], I[d], tmp);
				d--;
			}
			c--;
		}
		if (b > c) {
			break;
		}
		swap(I[b], I[c], tmp);
		b++;
		c--;
	}

	s = a - lowestPos;
	t = b - a;
	s = min(s, t);
	for (l = lowestPos, m = b - s; m < b; l++, m++) {
		swap(I[l], I[m], tmp);
	}

	s = d - c;
	t = highestPos - d;
	s = min(s, t);
	for (l = b, m = highestPos - s + 1; m <= highestPos; l++, m++) {
		swap(I[l], I[m], tmp);
	}

	s = b - a;
	t = d - c;
	if (s > 0) {
		QSufSortSortSplit(V, I, lowestPos, lowestPos + s - 1, numSortedChar);
	}

	// Update group number for equal portion
	a = lowestPos + s;
	b = highestPos - t;
	if (a == b) {
		// Sorted group
		V[I[a]] = a;
		I[a] = -1;
	} else {
		// Unsorted group
		for (c=a; c<=b; c++) {
			V[I[c]] = b;
		}
	}

	if (t > 0) {
		QSufSortSortSplit(V, I, highestPos - t + 1, highestPos, numSortedChar);
	}

}

/* Algorithm by Bentley & McIlroy.*/
static int QSufSortChoosePivot(int* __restrict V, int* __restrict I, const int lowestPos, 
							   const int highestPos, const int numSortedChar) {

	int m;
	int keyl, keym, keyn;
	int key1, key2, key3;
	int s;
	int numItem;

	numItem = highestPos - lowestPos + 1;

	m = lowestPos + numItem / 2;

	s = numItem / 8;
	key1 = KEY(V, I, lowestPos, numSortedChar);
	key2 = KEY(V, I, lowestPos+s, numSortedChar);
	key3 = KEY(V, I, lowestPos+2*s, numSortedChar);
	keyl = med3(key1, key2, key3);
	key1 = KEY(V, I, m-s, numSortedChar);
	key2 = KEY(V, I, m, numSortedChar);
	key3 = KEY(V, I, m+s, numSortedChar);
	keym = med3(key1, key2, key3);
	key1 = KEY(V, I, highestPos-2*s, numSortedChar);
	key2 = KEY(V, I, highestPos-s, numSortedChar);
	key3 = KEY(V, I, highestPos, numSortedChar);
	keyn = med3(key1, key2, key3);

	return med3(keyl, keym, keyn);


}

/* Quadratic sorting method to use for small subarrays. */
static void QSufSortInsertSortSplit(int* __restrict V, int* __restrict I, const int lowestPos, 
									const int highestPos, const int numSortedChar) {

	int i, j;
	int tmpKey, tmpPos;
	int numItem;
	int key[INSERT_SORT_NUM_ITEM], pos[INSERT_SORT_NUM_ITEM];
	int negativeSortedLength;
	int groupNum;

	numItem = highestPos - lowestPos + 1;

	for (i=0; i<numItem; i++) {
		pos[i] = I[lowestPos + i];
		key[i] = V[pos[i] + numSortedChar];
	}

	for (i=1; i<numItem; i++) {
		tmpKey = key[i];
		tmpPos = pos[i];
		for (j=i; j>0 && key[j-1] > tmpKey; j--) {
			key[j] = key[j-1];
			pos[j] = pos[j-1];
		}
		key[j] = tmpKey;
		pos[j] = tmpPos;
	}

	negativeSortedLength = -1;

	i = numItem - 1;
	groupNum = highestPos;
	while (i > 0) {
		I[i+lowestPos] = pos[i];
		V[I[i+lowestPos]] = groupNum;
		if (key[i-1] == key[i]) {
			negativeSortedLength = 0;
		} else {
			if (negativeSortedLength < 0) {
				I[i+lowestPos] = negativeSortedLength;
			}
			groupNum = i + lowestPos - 1;
			negativeSortedLength--;
		}
		i--;
	}

	I[lowestPos] = pos[0];
	V[I[lowestPos]] = groupNum;
	if (negativeSortedLength < 0) {
		I[lowestPos] = negativeSortedLength;
	}	

}
