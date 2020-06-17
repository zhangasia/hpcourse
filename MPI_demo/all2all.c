#include "myhead.h"
//用send,recv 等实现MPI_Alltoall
void all2all(sendbuf, sendcount, sendtype, recvbuf,recvcount, recvtype, comm, iam, np)
MPI_Comm comm;
int iam, np, sendcount, recvcount;
float* sendbuf, * recvbuf;
MPI_Datatype sendtype, recvtype;
{
	MPI_Status st;
	//int front, next;
	for (int i = 0; i < np; i++)
	{
		if (iam == i)
		{
			MPI_Sendrecv(&sendbuf[sendcount * i], sendcount, sendtype, i, 100,
				&recvbuf[recvcount * i], recvcount, recvtype, i, 100, comm, &st);
		}
		if (iam != i)
		{
			MPI_Send(&sendbuf[sendcount * i], sendcount, sendtype, i, i, comm);
			MPI_Recv(&recvbuf[recvcount * i], recvcount, recvtype, i, iam, comm, &st);
		}
	}
	return;
}
