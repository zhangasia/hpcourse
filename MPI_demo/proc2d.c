#include "myhead.h"
/*
generate 2-d processes produce row and column communicator for each 
process,according to the row major manner.The process array is p times q
*/
void proc2d(comm, np, iam, p, q, rowcom, colcom, rowid, colid)
//根据comm生成rowcom,colcom,根据iam生成rowid,colid
MPI_Comm comm, * rowcom, * colcom;
int np, iam, p, q, *rowid, *colid;
{
	int color, key, pxq;
	MPI_Comm valcom;
	pxq = p * q;
	if (np < pxq) return;//进程数不够

	if (iam < pxq) color = 0;
	else color = MPI_UNDEFINED;//通信子为空

	key = iam;
	MPI_Comm_split(comm, color, key, &valcom);

	if (valcom == MPI_COMM_NULL) return;

	//from row communicator
	color = iam / q;
	MPI_Comm_split(valcom, color, key, rowcom);

	// produce column communicator
	color = iam % q;
	MPI_Comm_split(valcom, color, key, colcom);

	MPI_Comm_rank(*colcom, rowid);
	MPI_Comm_rank(*rowcom, colid);

	return;
}