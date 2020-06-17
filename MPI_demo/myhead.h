#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
typedef struct {int a;float b[2];char c[3];}abc;
typedef struct { float a; int m; }floatint;
void mybegin(int*, char***, MPI_Comm*, int*, int*);
void myend();
void ring(int, int*, MPI_Comm, int, int);
//void datatype(char, MPI_Datatype*);
MPI_Datatype datatype(char);
//void mpistruct(MPI_Datatype*);
MPI_Datatype mpistruct();
//void diagonal(int, int, int, MPI_Datatype*, MPI_Datatype*);
MPI_Datatype diagonal(int, int, int);
//æÿ’Û≥À
/*
void init_a(int, int, int lda, float [][lda], int);
void init_b(int, int, int ldb, float [][ldb], int);
void matmul(int, int, int, int lda, float[][lda], \
	int ldb, float[][ldb], int ldc, float[][ldc]);
void rcmatmul(int, int, int, int, int, int, int lda, float[][lda], \
	int ldb, float[][ldb, int ldc, float[][ldc], int ldw, float[][ldw]);
*/
void init_a(int, int, int, float *, int);
void init_b(int, int, int, float *, int);
void matmul(int, int, int, int, float*, \
	int, float*, int, float*);
void rcmatmul(int, int, int, int, int, int, int, float*, \
	          int, float*, int, float*, int, float*);

void proc2d(MPI_Comm, int, int, int, int, MPI_Comm*, MPI_Comm*, int*, int*);

void snglscan(MPI_Comm, int, float, int, float*);
void gemmv(int, int, float*, int, float*, float*);
void iteration(MPI_Comm, int, int, int, int, float*, int, float*, float*, int);
void gemmm(int, int, int, float*, int, float*, int, float*, int);
void typemat(int, int, int, MPI_Datatype*);
void scopy(int, int, float*, int, float*, int);
void setinittab(int, int, int, int, int, int, float*, int, float*, int);
void cannon(MPI_Comm, MPI_Comm, int, int, int, int, int, int, float*, \
	int, float*, int, float*, int, float*, int, float*, int);

void all2all(float*, int, MPI_Datatype, float*, int, MPI_Datatype, MPI_Comm, int, int);

void inita(int, int, int,int*);
void initx(int, int*);
void gmv(int, int*, int*,int *);
void mv(MPI_Comm, int*, int*, int*,int*,int*, int, int, int, int);
void cpy(int, int*, int*);
void mysum(int, int, int, float*, int, float*);
