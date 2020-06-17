#include "myhead.h"
#include <string.h>

#define maxnp 16
void main(argc, argv)
int argc;
char** argv;
{
	MPI_Comm comm;
	int np, iam;
	int m, n,narray[100],marray[101];//narray,marray 验证自定义数据类型
	int matr[11][25];//MPI_Type_vector传输矩阵
	MPI_Datatype newtp ;
	MPI_Status st;
	abc x[10];
	MPI_Aint sizeabc, extnewtp,lowsize;
	float a[31][57];//验证diagonal
	float b[53][59], c[31][61], w[51][53],u[37][41];//矩阵乘
	float rhs[31], xx[31] = {0};//rhs 右端项
	int rcounts[maxnp];
	MPI_Aint displs[maxnp];
	FILE* fp;

	int i,j,en,k;
	int bnp[5];//从文件中读取的数据存放在该数组里

	int offside;

	floatint mxl, resmxl;

	int rowid, colid;
	MPI_Comm rowcom, colcom;

	MPI_Group grp1, newgrp;
	int ranks[10],gnp,giam;
	int p;

	int aaa[44][44], xxx[44], bbb[44],yyy[44],www[44];


#define zhangmv
	//start the MPI environment
	mybegin(&argc, &argv, &comm, &np, &iam);


#ifdef zhangmv
	m = 11;
	n = 44;
	if (iam == 0)
	{
		for (i = 0; i < m; i++)
		{
			if (i % 2 == 0) bbb[i] = 1;
			else bbb[i] = -1;
		}
	}
	MPI_Bcast(bbb, m, MPI_INT, 0, comm);
	inita(m, n,iam, aaa);
	initx(m, xxx);
	mv(comm, aaa, &xxx[m], bbb, yyy, www, 11, 44, np, iam);
	printf("\n x = %d,%d,%d on process %d\n", xxx[0], xxx[1], xxx[2], iam);
#endif
#ifdef zhangall2all
	for (int i = 0; i < 31; i++)
		for (int j = 0; j < 57; j++)
			a[i][j] = i + j;
//	MPI_Alltoall(a, 1, MPI_FLOAT, a, 1, MPI_FLOAT, comm);
	all2all(a, 1, MPI_FLOAT, b,1, MPI_FLOAT, comm, iam, np);
	printf("\n a = %f,%f,%f on process %d", b[0][0], b[0][1], b[0][2], iam);
#endif

#ifdef zhangcannon
	p = 3;
	if (np < 9) return;
	proc2d(comm, np, iam, p, p, &rowcom, &colcom, &rowid, &colid);
	/*
	if (iam == 0)
	{
		fp = fopen("inputmkn.txt", "r");
		i = fscanf(fp, "%*[^\n%*c");
		i = fscanf(fp, "%*[^\n%*c %d,%d,%d", &bnp[0], &bnp[1], &bnp[2]);
		fclose(fp);
		printf("nm = %d, k = %d, n = %d\n", bnp[0], bnp[1], bnp[2]);
	}
	MPI_Bcast(bnp, 3, MPI_INT, 0, comm);*/
	m = 11;
	k = 10;
	n = 12;
	if (iam < 9)
	{
		setinittab(p, rowid, colid, m, k, n, &a[0][0], 57, &b[0][0], 59);
	//	printf("\n a = %f,%f,%f,%f on process %d\n", \
			a[1][1], a[1][2], a[1][3], a[1][4], iam);
	//	printf("\n b = %f,%f,%f,%f on process %d\n", \
			b[1][1], b[1][2], b[1][3], b[1][4], iam);
		cannon(rowcom, colcom, p, rowid, colid, m, k, n, &a[0][0], 57, \
			& b[0][0], 59, &c[0][0], 61, &w[0][0], 53, &u[0][0], 41);
		printf("\n c = %f,%f,%f,%f on process %d\n", \
			c[1][1], c[1][2], c[1][3], c[1][4], iam);
	}

#endif

#ifdef zhanggroup
	MPI_Comm_group(comm, &grp1);
	ranks[0] = 1;
	ranks[1] = 3;
	MPI_Group_excl(grp1, 2, ranks, &newgrp);
	//MPI_Group_incl(grp1, 2, ranks, &newgrp);//把grap1里的1，3组成新的newgrp
	MPI_Group_size(newgrp, &gnp);
	MPI_Group_rank(newgrp, &giam);
	printf("\n The process %d,group %d\n", iam, giam);
#endif

#ifdef zhangiteration
	/*
	Ax = b,where A is diagonal,aii = 1/2,bi = i
	*/
	en = 5;
	n = en * np;
	for (i = 0; i < n; i++) rhs[i] = i;
	offside = iam * en;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < en; j++)
		{
			a[i][j] = 0.0;
			if (i == (j + offside)) a[i][j] = 0.5;
		}
	}
	/*
	Jacobi iteration x = (I-A)x+b
	*/
	iteration(comm, np, iam, n, en, &a[0][0], 57, rhs, xx, 10);
//	printf("\n x = %f,%f,on process %d\n", a[0][0], a[1][0], iam);
	printf("\n x = %f,%f,on process %d\n", xx[0], xx[1], iam);

#endif

#ifdef zhangsngl
	a[0][0] = iam + 1.0;
	snglscan(comm, iam, a[0][0], 2, &b[0][0]);
	if (iam == 2)
		printf("\nsum = %f\n", b[0][0]);
	MPI_Scan(&a[0][0], &b[0][0],1, MPI_FLOAT, MPI_SUM, comm);
	printf("\nEach process value = %f,in %d", b[0][0], iam);

#endif

#ifdef zhangproc2d
	if (np < 12) return;
	proc2d(comm, np, iam, 3, 4, &rowcom, &colcom, &rowid, &colid);
	printf("\nProcess %d = (%d,%d)\n", iam, rowid, colid);
#endif

#ifdef zhangreduce
	//MPI_Reduce-----------------------------------------------------------
	/*
	for (int i = 0; i < 31; i++)
		for (int j = 0; j < 57; j++)
			a[i][j] = i + j;
	MPI_Reduce(&a[1][1], &a[0][0], 1, MPI_FLOAT, MPI_SUM, 0, comm);
	*/
	
	//在Windows上，因为narray[]和marray[]是int类型，用%f打印输出为0
	//MPI_2INT是传递的两个数，和发送那一个的连着
	/*
	narray[0] = (iam + 1) * 20;
	narray[1] = iam;
//	MPI_Reduce(narray, marray, 1, MPI_2INT, MPI_MAXLOC, 0, comm);
	MPI_Allreduce(narray, marray, 1, MPI_2INT, MPI_MAXLOC, comm);
	printf("\nA = %d,location = %d in process %d \n", marray[0], marray[1], iam);
	*/
	//用结构体类型求最大值
	mxl.a = (iam + 1) * 20;
	mxl.m = iam;
	MPI_Reduce(&mxl, &resmxl, 1, MPI_FLOAT_INT, MPI_MAXLOC, 0, comm);
	if(iam == 0)
		printf("\nMax = %f,location = %d,in process %d", resmxl.a, resmxl.m, iam);



#endif

#ifdef zhangmatmul
	/*if (iam == 0)
	{
		fopen_s(&fp,"inputmkn.txt", "r");
		i = fscanf_s(&fp, "%*[^\n]%*c");
		i = fscanf_s(&fp, "%*[^\n]%*c %d,%d,%d", &bnp[0], &bnp[1], &bnp[2]);
		fclose(fp);
		printf("nm = %d, k = %d, n = %d\n", bnp[0], bnp[1], bnp[2]);
	}
	MPI_Bcast(bnp, 3, MPI_INT, 0, comm);
	m = bnp[0];
	k = bnp[1];
	n = bnp[2];*/
	m = 11;
	k = 45;
	n = 12;
	init_a(m, k, 57, a, iam);
	init_b(k, n, 59, b, iam);
//	printf("\n c = %f,%f,%f,%f on process %d\n", a[0][1], a[0][2], a[0][3], a[0][4], iam);
	rcmatmul(comm, np, iam, m, k, n, 57, a, 59, b, 61, c, 53, w);
	printf("\n c = %f,%f,%f,%f on process %d\n", c[1][1], c[1][2], c[1][3], c[1][4], iam);
#endif

#ifdef zhangbcast
	//MPI_Bcast--------------------------------------------------------------
	/*
	if (iam == 0)
	{
		for (int i = 0; i < 31; i++)
			for (int j = 0; j < 57; j++)
				a[i][j] = i + j;
		
	}
	MPI_Bcast(&a[0][0], 5, MPI_FLOAT, 0, comm);
	printf("\nA = %f,%f,%f in process %d\n", a[0][0], a[0][1], a[0][2], iam);
	*/
	//MPI_Gather-------------------------------------------------------------
	/*
	for (int i = 0; i < 31; i++)
		for (int j = 0; j < 57; j++)
			a[i][j] = i + j;
	int j = iam * 5;
	MPI_Gather(&a[0][j], 3, MPI_FLOAT, &a[1][0], 3, MPI_FLOAT,0, comm);

	printf("\nA = %f,%f,%f,%f,%f,%f in process %d\n",\
		a[1][0], a[1][1], a[1][2], a[1][3], a[1][4], a[1][5], iam);
	*/
	//MPI_Gatherv-----------------------------------------------------------
	/*
	for (int i = 0; i < 31; i++)
		for (int j = 0; j < 57; j++)
		{
			if (iam == 0)
			{
				a[i][j] = 0;
			}
			else
			{
				a[i][j] = i + j;
			}	
		}
			
	for (int i = 0; i < np; i++)
	{
		rcounts[i] = 3;
		displs[i] =  i * 5;//每隔5个放一个
	}
	int j = iam * 5;
	MPI_Gatherv(&a[0][j], 3, MPI_FLOAT, &a[1][0], rcounts, displs, MPI_FLOAT, 0, comm);
	if (iam == 0)
	{
		printf("\nA = %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f in process %d\n", \
			a[1][0], a[1][1], a[1][2], a[1][5], a[1][6], a[1][7],\
			a[1][10], a[1][11], a[1][12], a[1][15], a[1][16], a[1][17],
			iam);
	}
	*/
	//MPI_Scatter----------------------------------------------------------
	/*
	for (int i = 0; i < 31; i++)
		for (int j = 0; j < 57; j++)
			a[i][j] = i + j;
	int j = iam * 5;
	MPI_Gather(&a[0][j], 3, MPI_FLOAT, &a[1][0], 3, MPI_FLOAT, 0, comm);
	printf("\nA = %f,%f,%f,%f,%f,%f in process %d\n", \
		a[1][0], a[1][1], a[1][2], a[1][3], a[1][4], a[1][5], iam);
	MPI_Scatter(&a[1][0], 3, MPI_FLOAT, &a[0][0], 3, MPI_FLOAT, 0, comm);
	printf("\nA = %f,%f,%f in process %d\n",a[0][0], a[0][1], a[0][2],iam);
	*/
	//MPI_Allgather-----------------------------------------------------------------
	/*
	if (iam == 0)
	{
		for (int i = 0; i < 31; i++)
			for (int j = 0; j < 57; j++)
				a[i][j] = i + j;

	}
	MPI_Bcast(&a[0][0], 5 * np, MPI_FLOAT, 0, comm);
//	printf("\nA = %f,%f,%f in process %d\n", a[0][0], a[0][1], a[0][2], iam);
	MPI_Scatter(&a[0][0], 5, MPI_FLOAT, &a[1][0], 5, MPI_FLOAT, 1, comm);
//	printf("\nAA = %f,%f,%f,%f,%f in process %d\n", a[1][0], a[1][1], a[1][2], a[1][3], a[1][4], iam);
	MPI_Allgather(&a[1][0], 5, MPI_FLOAT, &a[2][0], 5, MPI_FLOAT, comm);
//	printf("\nA = %f,%f,%f,in process %d\n", a[2][6], a[2][1], a[2][10], iam);
//	printf("\nAA = %f,%f,%f,%f,%f in process %d\n", a[1][0], a[1][1], a[1][2], a[1][3], a[1][4], iam);
	MPI_Alltoall(&a[0][0], 5, MPI_FLOAT, &a[3][0], 5, MPI_FLOAT, comm);
	printf("\nA = %f,%f,%f,%f,%f,in process %d\n", a[0][0], a[0][1], a[0][2], a[0][3], a[0][4], iam);
	printf("\nA = %f,%f,%f,%f,%f in process %d\n", a[3][0], a[3][1], a[3][2], a[3][3], a[3][4], iam);
	*/
	for (int i = 0; i < 31; i++)
		for (int j = 0; j < 57; j++)
			a[i][j] = iam*10;
	MPI_Alltoall(&a[0][0], 1, MPI_FLOAT, &a[3][0], 1, MPI_FLOAT, comm);
  //	printf("\nAA = %f,%f,%f,%f,%f,in process %d\n", a[0][0], a[0][1], a[0][2], a[0][3], a[0][4], iam);
	printf("\nA = %f,%f,%f,%f,%f in process %d\n", a[3][0], a[3][1], a[3][2], a[3][3], a[3][4],  iam);

#endif

#ifdef zhangdiagonal
	newtp = diagonal(2, 3, 57);
	MPI_Type_get_extent(newtp, &lowsize, &extnewtp);
	printf("\n lowbound = %ld extent=%ld\n", lowsize, extnewtp);//57*2 * 4
	if (iam == 0)
	{
		for (int i = 0; i < 31; i++)
			for (int j = 0; j < 57; j++)
				a[i][j] = i + j;
		MPI_Send(a, 3, newtp, 1, 5, comm);
	}
	else if (iam == 1)
	{
		MPI_Recv(a, 3, newtp, 0, 5, comm, &st);
		printf("\nA=%f,%f,%f,%f,%f,%f,%f,%f\n", \
			a[0][0], a[0][1], a[1][0], a[1][1], \
			a[2][0], a[2][1], a[2][3], a[2][4]);
	}

	MPI_Type_free(&newtp);for (int i = 0; i < 31; i++)
		for (int j = 0; j < 57; j++)
			a[i][j] = i + j;
#endif

#ifdef zhangstruct
	newtp = mpistruct();
	if (iam == 0)
	{
		for (m = 0; m < 10; m++)
		{
			x[m].a = m;
			x[m].b[0] = 20.0 * (m + 1);
			x[m].b[1] = 30.0 * (m + 1);
			x[m].c[0] = 'a' + 3 * m;
			x[m].c[1] = 'b' + 3 * m;
			x[m].c[2] = 'c' + 3 * m;
		}
		MPI_Send(x, 3, newtp, 1, 5, comm);
		sizeabc = sizeof(abc);
		MPI_Type_get_extent(newtp, &lowsize,&extnewtp);
		printf("\n sizeof=%ld, and lowbound = %ld extent=%ld\n", sizeabc,lowsize ,extnewtp);
	}
	if (iam == 1)
	{
		MPI_Recv(x, 3, newtp, 0, 5, comm, &st);
		printf("\n values are %d,%f,%f,%c,%c,%c\n", \
			x[0].a, x[0].b[0], x[0].b[1], x[0].c[0], x[0].c[1], x[0].c[2]);
		printf("\n values are %d,%f,%f,%c,%c,%c\n", \
			x[1].a, x[1].b[0], x[1].b[1], x[1].c[0], x[1].c[1], x[1].c[2]);
	}

#endif

#ifdef zhangdatatype
	for (m = 0; m < 100; m++)
	{
		narray[m] = m;
	//	marray[m] = 1;
	}
	//MPI_Type_contiguous
	/*
	newtp = datatype('c');
	if (iam == 0)
	{
		MPI_Send(narray, 3, newtp, 1, 5, comm);
	}
	if (iam == 1)
	{
	//	MPI_Recv(marray, 3, newtp, 0, 5, comm, &st);
		MPI_Recv(marray, 6, MPI_INT, 0, 5, comm, &st);//也可以直接接收6个MPI_INT
		printf("\nData on Process %d are %d,%d,%d,%d,%d,%d\n", \
			iam, marray[0], marray[1], marray[2], marray[3], marray[4], marray[5]);
	}
	*/
	//MPI_Type_vector
	
	newtp = datatype('v');
	if (iam == 0)
	{
		for (m = 0; m < 10; m++)
		{
			for (n = 0; n < 25; n++)
			{
				matr[m][n] = m + n;
			}
		}
		MPI_Send(narray, 1, newtp, 1, 5, comm);
		//MPI_Send(matr, 1, newtp, 1, 5, comm);
	}
	if (iam == 1)
	{

		MPI_Recv(marray, 1, newtp, 0, 5, comm, &st);
		printf("\nData on Process %d are %d,%d,%d,%d,%d,%d\n", \
			iam, marray[0], marray[1], marray[2], marray[3], marray[6], marray[7]);
		//MPI_Recv(matr, 1, newtp, 0, 5, comm, &st);
		//printf("\nData on Process %d are %d,%d,%d,%d,%d,%d\n", \
			iam, matr[0][0], matr[0][1], matr[0][2],\
				matr[1][0], matr[1][1], matr[1][2]);
	}
	
	//MPI_Type_indexed
	/*
	newtp = datatype('i');
	if (iam == 0)
	{
		MPI_Send(narray, 1, newtp, 1, 5, comm);
	}
	if (iam == 1)
	{
		MPI_Recv(marray, 3, newtp, 0, 5, comm, &st);
		
		printf("\nData on Process %d are %d,%d,%d,%d,%d\n", \
			iam, marray[0], marray[1], marray[2], marray[5], marray[6]);
	}
	*/
		
	//printf("\nData type is created\n");
	MPI_Type_free(&newtp); matr[0][0],

#endif

#ifdef zhangring
	//main body here
	m = iam;
	n = 100;
	ring(m, &n, comm, np, iam);
	printf("\nIn process %d n = %d!\n", iam,n);
#endif

	myend();
}