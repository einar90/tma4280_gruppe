#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "solvers.h"
#include "poissoncommon.h"

double alpha=0.0;
double tol=1e-6;

double exact(double x, double y)
{
  return x*(pow(x,5)-1.0)*y*(pow(y,5)-1.0);
}

double source(double x, double y)
{
  return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
}





void DiagonalizationPoisson2Dfst(Matrix b, const Vector lambda)
{
  int i,j;
  Matrix ut = cloneMatrix(b);
  Vector buf = createVector(4*(b->rows+1));
  int N=b->rows+1;
  int NN=4*N;

#pragma omp parallel for schedule(static) private(i)
  for (i=0;i<b->cols;++i)
    fst(b->data[i], &N, buf->data, &NN);

  transposeMatrix(ut, b);

#pragma omp parallel for schedule(static) private(i)
  for (i=0;i<ut->cols;++i)
    fstinv(ut->data[i], &N, buf->data, &NN);

  for (j=0;j<b->cols;++j)
    for (i=0;i<b->rows;++i)
      ut->data[j][i] /= (lambda->data[i]+lambda->data[j]+alpha);

#pragma omp parallel for schedule(static) private(i)
  for (i=0;i<b->cols;++i)
    fst(ut->data[i], &N, buf->data, &NN);
  transposeMatrix(b, ut);

#pragma omp parallel for schedule(static) private(i)
  for (i=0;i<ut->cols;++i)
    fstinv(b->data[i], &N, buf->data, &NN);

  freeMatrix(ut);
  freeVector(buf);
}






int main(int argc, char** argv)
{
  int i, j, N, flag, local;
  Matrix b, e;
  Vector grid, lambda=NULL;
  double time, sum, h;

  flag=12;
  local = 0;

  N=atoi(argv[1]);
  if (N < 0) {
    printf("invalid problem size given\n");
    return 2;
  }
  if ((N & (N-1)) != 0) {
    printf("need N to be a power-of-two for fst-based diagonalization\n");
    return 5;
  }

  // Generate 1D grid from 0 to 1
  grid = equidistantMesh(0.0, 1.0, N);

  // Create matrices
  b = createMatrix(N-1,N-1);
  e = createMatrix(N-1,N-1);

  // Populate b with the source function evaluated for grid values
  evalMeshInternal2(b, grid, source, local);

  h = grid->data[1]-grid->data[0];
  scaleVector(b->as_vec, pow(h, 2));
  evalMeshInternal2(e, grid, exact, local);
  axpy(b->as_vec, e->as_vec, alpha);

  // Generate eigenvalues
  lambda = generateEigenValuesP1D(N-1);

  time = WallTime();

  // Solving system b
  DiagonalizationPoisson2Dfst(b, lambda);

  printf("elapsed: %f\n", WallTime()-time);


  axpy(b->as_vec,e->as_vec,-1.0);

  // // Printing b
  // for (i = 0; i < N; ++i)
  // {
  //   for (j = 0; j < N; ++j)
  //   {
  //     printf("%4.2f ", b->data[i][j]);
  //   }
  //   printf("\n");
  // }

  printf("max error: %e\n", maxNorm(b->as_vec));



  freeMatrix(b);
  freeMatrix(e);
  freeVector(grid);
  if (lambda)
    freeVector(lambda);
  return 0;
}
