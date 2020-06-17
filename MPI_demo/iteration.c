#include "myhead.h"

void iteration(comm, np, iam, n, en, a, lda, b, x, num)
//ÿ�������Ͼ����������n/np.en ����ÿ�������е�����,num �����������
MPI_Comm comm;
int np, iam, n, en, lda, num;
float* a, * b, * x;
{
	int i, j, *rc;
	float *y;
	rc = (int* )malloc(np * sizeof(int));
	for (i = 0; i < np; i++) rc[i] = en;//���տ�Ĵ�С
	y = (float* )malloc(n * sizeof(float));//y = ax + b;
	for (i = 0; i < num; i++)
	{
		if (iam == 0) //����iam = 0��ʱ����b,�������û��b
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