#include "myhead.h"

/*
matrix multiplication by a vector
*/

//ax + y = y
void gemmv(m, n, a, lda, x, y)
int m, n, lda;
float *a, *x, *y;
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			y[i] += a[i * lda + j] * x[j];
	return;
}