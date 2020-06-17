#include "myhead.h"
//行列分块矩阵乘
/*
The matrix A is partitioned by row and B by column,a_ij = i + j
b_ij = 1 if J is even, else -1
*/

// m X k is the block matrix order in iam,the full matrix A is np * m X k mareix
//把A矩阵分成了np个m*k的小矩阵

void init_a(m, k, lda, a, iam)
int m, k, lda, iam;
float *a;
{
	int i, j, offside;
	offside = iam * m;
	for (i = 0; i < m; i++)
		for (j = 0; j < k; j++)
			a[i*lda + j] = offside + i + j;
	return;
}
// k X n is the order of matrix B in iam,original B is k * np * n
void init_b(k, n, ldb, b, iam)
int k, n, ldb, iam;
float *b;
{
	int i, j ,offside;
	offside = iam * iam;
	for (i = 0; i < k; i++)
		for (j = 0; j < n; j++)
			b[i * ldb + j] = 1.0 - 2.0 * ((i + j + offside) % 2);
	return;
}
//c = a X b
/*
void matmul(m, k, n, lda, a, ldb, b, ldc, c)
int m, k, n, lda, ldb, ldc;
float a[][lda], b[][ldb], c[][ldc];
{
	int i, j, l;
	for(i = 0;i < m;i++)
		for (j = 0; j < n; j++)
		{
			c[i][j] = 0.0;
			for (l = 0; l < k; l++)
				c[i][j] += a[i][l] + b[l][j];
		}
	return;
}
*/
void matmul(m, k, n, lda, a, ldb, b, ldc, c)
int m, k, n, lda, ldb, ldc;
float *a, *b, *c;
{
	int i, j, l;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			c[i*ldc + j] = 0.0;
			for (l = 0; l < k; l++)
				c[i * ldc + j] += a[i * lda + l] * b[l * ldb + j];
		}
	return;
}

//row-column partitioned algorithm
void rcmatmul(comm, np, iam, m, k, n, lda, a, ldb, b, ldc, c, ldw, w)
//w 临时空间，send b,recv w。m是原矩阵的行数除以np,n(原矩阵n * np)列数除以np.k为a列数，b行数.
MPI_Comm comm;
int np, iam, m, k, n, lda, ldb, ldc, ldw;
//float a[][lda], b[][ldb], c[][ldc], w[][ldw];
float *a, *b, *c, *w;
{
	int i, front, next, l;
	MPI_Datatype rectb, rectw;
	MPI_Status st;
	//create a new datatype for matrix b
	MPI_Type_vector(k, n, ldb, MPI_FLOAT, &rectb);
	MPI_Type_vector(k, n, ldw, MPI_FLOAT, &rectw);
	MPI_Type_commit(&rectb);
	MPI_Type_commit(&rectw);
	l = iam * n;
	front = (np + iam - 1) % np;
	next = (iam + 1) % np;
	for (i = 0; i < np - 1; i++)
	{
		if (i % 2 == 0)
		{
			matmul(m, k, n, lda, a, ldb, b, ldc, &c[l]);
			MPI_Sendrecv(b, 1, rectb, front, 1, w, 1, rectw, next, 1, comm, &st);
		}
		else
		{
			matmul(m, k, n, lda, a, ldw, w, ldc, &c[l]);
			MPI_Sendrecv(w, 1, rectw, front, 1, b, 1, rectb, next, 1, comm, &st);
		}
		l += n;
		if (l == np * n) l = 0;
	}
	if ((np - 1) % 2 == 0)
		matmul(m, k, n, lda, a, ldb, b, ldc, &c[l]);
	else
		matmul(m, k, n, lda, a, ldw, w, ldc, &c[l]);
	MPI_Type_free(&rectb);
	MPI_Type_free(&rectw);
	return;
}