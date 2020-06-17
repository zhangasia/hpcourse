#include "myhead.h"
void ring(m, n, comm, np, iam)
int m, * n, np, iam;
MPI_Comm comm;
{
	int front, next;
	MPI_Status st,sts[2];
	MPI_Request reqs[2], sreq, rreq;
	front = (np + iam - 1) % np; //防止出现负数
	next = (iam + 1) % np;
	/*
	if (iam % 2 == 0)
	{
		MPI_Send(&m, 1, MPI_INT, next, 1, comm);
		MPI_Recv(n, 1, MPI_INT, front, 1, comm, &st);

	}
	else
	{
		MPI_Recv(n, 1, MPI_INT, front, 1, comm, &st);
		MPI_Send(&m, 1, MPI_INT, next, 1, comm);
	}
	*/
	// 可能会造成死锁
	/*
	if (iam == 0)
	{
		MPI_Send(&m, 1, MPI_INT, next, 1, comm);
	}
	else if (iam = n - 1)
	{
		MPI_Recv(n, 1, MPI_INT, front, 1, comm, &st);
	}
	else
	{
		MPI_Sendrecv(&m, 1, MPI_INT, next, 1, n, 1, MPI_INT, front, 1, comm, &st);
	}
	*/
	//0进程不接收，np-1进程不发
	/*
	if (iam == 0) front = MPI_PROC_NULL;
	if (iam == np - 1) next = MPI_PROC_NULL;
	MPI_Sendrecv(&m, 1, MPI_INT, next, 1, n, 1, MPI_INT, front, 1, comm, &st);
	*/
	//测试非阻塞send和recv,用waitall
	
	MPI_Isend(&m, 1, MPI_INT, next, 1, comm, &reqs[0]);
	MPI_Irecv(n, 1, MPI_INT, front, 1, comm, &reqs[1]);
	MPI_Waitall(2, reqs, sts);
	//MPI_Request_free(reqs);//报错
	
	//用分开wait
	/*
	MPI_Isend(&m, 1, MPI_INT, next, 1, comm, &sreq);
	MPI_Irecv(n, 1, MPI_INT, front, 1, comm, &rreq);
	MPI_Wait(&sreq, &st);
	MPI_Wait(&rreq, &st);
	*/
	return;
}