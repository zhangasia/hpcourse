#include "myhead.h"
//c = c + a*b(¾ØÕó)
void gemmm(m, k, n, a, lda, b, ldb, c, ldc)
int m, k, n;//¾ØÕóÎ¬Êı
int lda, ldb, ldc;
float * a, * b, * c;
{
	int i, j, l;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			for (l = 0; l < k; l++)
				c[i * ldc + j] += a[i * lda + l] * b[l * ldb + j];
	return;
}
//define a matrix data type
void typemat(m, n, lda, newtp)
int m, n, lda;
MPI_Datatype* newtp;
{
	MPI_Type_vector(m, n, lda, MPI_FLOAT, newtp);
	return;
}
//copy a to b
void scopy(m, n, a, lda, b, ldb)
int m, n, lda, ldb;
float* a, * b;
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			b[i * ldb + j] = a[i * lda + j];
		}
	}
	return;
}
// A is a matrix,where aij = i + j,bij = 1,if (i + j) is even,else -1
void setinittab(p, myrow, mycol, m, k, n, a, lda, b, ldb)
int p, myrow, mycol, m, k, n, lda, ldb;
float* a, * b;
{
	int i, j, offsizea, offsizeb;
	offsizea = m * myrow + k * mycol;
	offsizeb = k * myrow + n * mycol;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			a[i * lda + j] = i + j + offsizea;
		}
	}
	for (i = 0; i < k; i++)
	{
		for (j = 0; j < n; j++)
		{
			/*
			if ((i + j + offsizeb) % 2 == 0)
			{
				b[i * ldb + j] = 1.0;
			}
			else
			{
				b[i * ldb + j] = -1.0;
			}*/
			b[i * ldb + j] = 1.0 - 2.0 * ((i + j + offsizeb) % 2);
		}
	}
}