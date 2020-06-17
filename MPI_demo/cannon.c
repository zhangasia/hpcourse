#include "myhead.h"
//Cannon Algorithm implemetation,the each A as m*k,and B is k*n,
//so the total matrix size is np*m*nq*k for matrix A,and forth.
//a[i][j] = a[i * lda + j]
void cannon(rowcom, colcom, p, myrow, mycol, m, k, n, a, lda, b, ldb, c, ldc, at, ldaw, bt, ldbw)
MPI_Comm rowcom, colcom;
int p, myrow, mycol, m, k, n, lda, ldb, ldc, ldaw, ldbw;
float* a, * b, * c, * at, * bt;//������һά����
{
	int i, j, l, front, next;
	MPI_Status st;
	MPI_Datatype btp, attp, bttp;//�����µ���������
   // typemat(m,k,lda,&atp);
   // MPI_Type_commit(&atp);
	typemat(k, n, ldb, &btp);
	MPI_Type_commit(&btp);
	typemat(m, k, ldaw, &attp);
	MPI_Type_commit(&attp);
	typemat(k, n, ldbw, &bttp);
	MPI_Type_commit(&bttp);

	l = myrow;
	front = (myrow - 1 + p) % p;
	next = (myrow + 1) % p;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			c[i * ldc + j] = 0.0;
		}
	}

	for (i = 0; i < p; i++)
	{
		if (mycol == l) //�Խ����ϵ�Ԫ��
		{
			scopy(m, k, a, lda, at, ldaw);//��a���Ƹ�at
		}
		MPI_Bcast(at, 1, attp, l, rowcom);//�Խ����ϵ�Ԫ�ع㲥���е�ÿ��Ԫ��
		gemmm(m, k, n, at, ldaw, b, ldb, c, ldc);
		if (i == p - 1) continue;
		//��b�����ƶ�
		MPI_Sendrecv(b, 1, btp, front, 1, bt, 1, bttp, next, 1, colcom, &st);
		scopy(k, n, bt, ldbw, b, ldb);
		l = (l + 1) % p;
	}
	return;
}