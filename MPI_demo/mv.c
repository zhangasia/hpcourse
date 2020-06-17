#include "myhead.h"
void inita(m, n,iam, a)
int m, n,iam;
int* a;
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			a[i * m + j] = i * n + j + iam * m;
		}
			
}
void initx(m, x)
int m;
int* x;
{
	int i;
	for (i = 0; i < m; i++)
	{
		if (i % 2 == 0) x[i] = 1;
		else x[i] = -1;
	}
}

void gmv(m, a, x, y)
int m;
int* a, * x, *y;
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
			y[i] += a[i * m + j] * x[j];
	return;
}
void cpy(m, x, y)
int m;
int* x, * y;
{
	for (int i = 0; i < m; i++)
		x[i] = y[i];
}
void mv(comm, a, x, b, y, w, m, n, np, iam)
MPI_Comm comm;
int m, n, np, iam;
int* a, * b, * x, *y,*w;
{
	int i, j, l;
	int front, next;
	MPI_Status st;
	front = (np + iam - 1) % np;
	next = (iam + 1) % np;
	l = 0;
	for (i = 0; i < m; i++)
	{
//		b[i] = x[i];
		y[i] = b[i];
		w[i] = x[i];

	}
	for (i = 0; i < np - 1; i++)
	{
		if (i % 2 == 0)
		{
			gmv(m, &a[l], x, y);
			MPI_Sendrecv(x, m, MPI_INT, front, 1, w, m, MPI_INT, next, 1, comm, &st);
		}
		else
		{
			gmv(m, &a[l], w, y);
			MPI_Sendrecv(w, m, MPI_INT, front, 1, x, m, MPI_INT, next, 1, comm, &st);
		}
		l += m * m;
		if (l == n * m) l == 0;
	}
	if ((np - 1) % 2 == 0)
	{
		gmv(m, &a[l], x, y);
	}
	else
	{
		gmv(m, &a[l], w, y);
	}
	cpy(m, x, y);
}