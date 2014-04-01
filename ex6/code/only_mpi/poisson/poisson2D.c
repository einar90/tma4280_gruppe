#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "poissoncommon.h"

double alpha=0.0;
double tol=1e-6;

double exact(double x, double y)
{
  return x*(pow(x,5)-1.0)*y*(pow(y,5)-1.0);
}

double source(double x, double y)
{
  // return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
  return 2*x+y;
}


Matrix createPartMatrix(int colsToSend, int N, Matrix b, int i, int* displacements)
{
  int c, r;
  Matrix copyMatrix = createMatrix(N-1, colsToSend);
  for(c = 0; c < colsToSend; c++) {
    for (r = 0; r < b->rows; ++r)
    {
      copyMatrix->data[c][r] = b->data[c + displacements[i-1]][r];
    }
  }
  return copyMatrix;
}


void distributeCols(int* partlens, int* displacements, int N, Matrix b,
                    Matrix copyMatrix, int processor_count, int tag)
{
  int i;
  for(i = 1; i < (processor_count); i++) {
    int colsToSend = partlens[i-1];
    copyMatrix = createPartMatrix(colsToSend, N, b, i, displacements);
    MPI_Send(&copyMatrix->as_vec->data[0], copyMatrix->as_vec->len, MPI_DOUBLE,
             i, tag, MPI_COMM_WORLD);
    freeMatrix(copyMatrix);
  }
}


void recvParts(int* partlens, int* displacements, Matrix b, int N, int processor_count, int tag)
{
  int i;
  for (i = 1; i < processor_count; ++i)
  {
    int first_index = displacements[i-1] * (N-1);
    int part_len = partlens[i-1] * (N-1);
    MPI_Recv(&b->as_vec->data[first_index], part_len, MPI_DOUBLE, i,
             tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}



void DiagonalizationPoisson2Dfst(Matrix b, int size, int rank,
                                 MPI_Status status, int *displacements,
                                 int *partlens)
{
  int i,j,k,r,c;
  int N=b->cols+1;
  Vector lambda;
  Matrix ut = cloneMatrix(b);
  Vector buf = createVector(4*(b->rows+1));
  int NN=4*N;
  double time;

  Matrix copyMatrix;

  if(rank == 0) {
    lambda = generateEigenValuesP1D(N-1);

    distributeCols(partlens, displacements, N, b, copyMatrix, size, 100);
  }

  if(rank != 0) {
    MPI_Recv(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE,
             0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    time = WallTime();
    for (i=0;i<b->cols;++i)
      fst(b->data[i], &N, buf->data, &NN);
    printf("Time spent on first fst: %f\n", WallTime() - time);

    MPI_Send(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE, 0, 200, MPI_COMM_WORLD);

    printf("b after recv for rank %d\n", rank);
    for (c = 0; c < b->cols; ++c)
    {
      for (r = 0; r < b->rows; ++r)
      {
        printf("%7.4f\t", b->data[c][r]);
      }
      printf("\n");
    }
  }


  if (rank == 0)
  {
    recvParts(partlens, displacements, b, N, size, 200);
    time = WallTime();
    transposeMatrix(ut, b);
    printf("Time spent on first transpose: %f\n", WallTime() - time);

    distributeCols(partlens, displacements, N, b, copyMatrix, size, 300);
  }


  if (rank != 0)
  {
    MPI_Recv(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE,
             0, 300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    time = WallTime();
    for (i=0;i<ut->cols;++i)
      fstinv(ut->data[i], &N, buf->data, &NN);
    printf("Time spent on first fstinv: %f\n", WallTime() - time);

    MPI_Send(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE, 0, 400, MPI_COMM_WORLD);
  }


  if (rank == 0)
  {
    recvParts(partlens, displacements, b, N, size, 400);
    time = WallTime();
    for (j=0;j<b->cols;++j){
      for (i=0;i<b->rows;++i){
        ut->data[j][i] /= (lambda->data[i]+lambda->data[j]+alpha);
      }
    }
    printf("Time spent computing lambdas: %f\n", WallTime() - time);
  }


  if(rank == 0) {

    time = WallTime();
    for (i=0;i<b->cols;++i)
      fst(ut->data[i], &N, buf->data, &NN);
    printf("Time spent on second fst: %f\n", WallTime() - time);

    time = WallTime();
    transposeMatrix(b, ut);
    printf("Time spent on second transpose: %f\n", WallTime() - time);  // freeMatrix(ut);

    time = WallTime();
    for (i=0;i<ut->cols;++i)
      fstinv(b->data[i], &N, buf->data, &NN);
    printf("Time spent on second fstinv: %f\n", WallTime() - time);

  }

  freeMatrix(ut);
  freeVector(buf);
  printf("Freed ut and buf on rank %d.\n", rank);
}





int main(int argc, char** argv)
{
  int i, j, N, flag, local, r, c;
  Matrix b, e;
  Vector grid;
  double time, sum, h;
  flag = 12;
  N=atoi(argv[1]);

  // MPI Vars
  int rank, size, tag;
  MPI_Status status;
  init_app(argc, argv, &rank, &size);

  int *partlens;
  int *displacements;
  splitVector(N-1, size-1, &partlens, &displacements);

  if (flag == 12 && (N & (N-1)) != 0) {
    printf("need N to be a power-of-two for fst-based diagonalization\n");
    return 5;
  }

  local = 0;

  // Master har ansvar for dette


  // Find a split size, exclude one processor as its master and handling send/recv

  if(rank == 0) {
    grid = equidistantMesh(0.0, 1.0, N);
    b = createMatrix(N-1,N-1);
    e = createMatrix(N-1,N-1);

    evalMeshInternal2(b, grid, source, local);
    h = grid->data[1]-grid->data[0];
    scaleVector(b->as_vec, pow(h, 2));
    evalMeshInternal2(e, grid, exact, local);
    axpy(b->as_vec, e->as_vec, alpha);

    time = WallTime();
  }

  if (rank != 0) {
    b = createMatrix(N-1, partlens[rank-1]);
  }



  DiagonalizationPoisson2Dfst(b, size, rank, status, displacements, partlens);
  printf("Diagonalization complete on rank %d\n", rank);


  if (rank == 0)
  {
    printf("elapsed: %f\n", WallTime()-time);
  }

  //axpy(b->as_vec,e->as_vec,-1.0);

  //printf("max error: %e\n", maxNorm(b->as_vec));

  // for (i = 0; i < b->rows; ++i)
  // {
  //   for (j = 0; j < b->cols; ++j)
  //   {
  //     printf("%f ", b->data[i][j]);
  //   }
  //   printf("\n");
  // }


  if(rank == 0) {
    freeMatrix(e);
    freeVector(grid);

    printf("b at end for rank %d\n", rank);
    for (r = 0; r < b->rows; ++r)
    {
      for (c = 0; c < b->cols; ++c)
      {
        printf("%7.4f\t", b->data[c][r]);
      }
      printf("\n");
    }
  }
  freeMatrix(b);

  close_app();
  return 0;
}
