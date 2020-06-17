#include "myhead.h"
MPI_Datatype diagonal(m, n, lda)
int m, n, lda;//获取小矩阵的m,n 为行列数,lda 为步长

{
	MPI_Datatype newtp,  rect;
	MPI_Datatype oldtype[2];
	int blken[2];
	MPI_Aint displs[2];
	MPI_Type_vector(m, n, lda, MPI_FLOAT, &rect);
	MPI_Type_commit(&rect);
	
	oldtype[0] = rect;
	oldtype[1] = MPI_UB;
	blken[0] = 1;
	blken[1] = 1;
	displs[0] = 0;
	displs[1] = sizeof(float) * (m * lda + n);
	MPI_Type_create_struct(2, blken, displs, oldtype, &newtp);
	MPI_Type_commit(&newtp);
	MPI_Type_free(&rect);
	return newtp;
}