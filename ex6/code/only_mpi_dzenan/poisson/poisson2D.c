#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "poissoncommon.h"




double source(double x, double y)
{
  return 5*x+y;
}


int main(int argc, char** argv)
{
  int rank, size;
  init_app(argc, argv, &rank, &size);
  MPI_Status status;

  int i, j, N, flag, local;
  Matrix b, e;
  Vector grid, lambda=NULL;
  double time, sum, h;
  flag = 12;
  N=atoi(argv[1]);

  if (flag == 12 && (N & (N-1)) != 0) {
    printf("need N to be a power-of-two for fst-based diagonalization\n");
    return 5;
  }

  local = 0;

  grid = equidistantMesh(0.0, 1.0, N);
  b = createMatrix(N-1,N-1);
  e = createMatrix(N-1,N-1);

  evalMeshInternal2(b, grid, source, local);


  printf("Before the matrix transpose\n");
  for (i = 0; i < b->rows; ++i)
  {
    for (j = 0; j < b->cols; ++j)
    {
      printf("%f ", b->data[i][j]);
    }
    printf("\n");
  }

  // Do Transpose
  if(rank == 0){
    printf("I am the master\n");
  } 
  else {
    printf("I am a slave\n");
  }











  printf("After matrix transpose\n");
  for (i = 0; i < b->rows; ++i)
  {
    for (j = 0; j < b->cols; ++j)
    {
      printf("%f ", b->data[i][j]);
    }
    printf("\n");
  }


  freeMatrix(b);
  freeMatrix(e);
  freeVector(grid);

  close_app();
  return 0;
}
