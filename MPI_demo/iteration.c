#include "myhead.h"

void iteration(comm, np, iam, n, en, a, lda, b, x, num)
//每个进程上矩阵的列数，n/np.en 代表每个进程中的列数,num 代表迭代次数
MPI_Comm comm;
int np, iam, n, en, lda, num;
float* a, * b, * x;
{
	int i, j, *rc;
	float *y;
	rc = (int* )malloc(np * sizeof(int));
	for (i = 0; i < np; i++) rc[i] = en;//接收块的大小
	y = (float* )malloc(n * sizeof(float));//y = ax + b;
	for (i = 0; i < num; i++)
	{
		if (iam == 0) //假设iam = 0的时候有b,其他情况没有b
			for (j = 0; j < n; j++) y[j] = b[j];
		else
			for (j = 0; j < n; j++) y[j] = 0.0;
		gemmv(n, en, a, lda, x, y);
		MPI_Reduce_scatter(y, x, rc, MPI_FLOAT, MPI_SUM, comm);
	}
	free(y);
	free(rc);
	return;
}