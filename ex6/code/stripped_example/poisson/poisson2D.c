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

  for (i=0;i<b->cols;++i)
    fst(b->data[i], &N, buf->data, &NN);
  transposeMatrix(ut, b);
  for (i=0;i<ut->cols;++i)
    fstinv(ut->data[i], &N, buf->data, &NN);

  for (j=0;j<b->cols;++j)
    for (i=0;i<b->rows;++i)
      ut->data[j][i] /= (lambda->data[i]+lambda->data[j]+alpha);

  for (i=0;i<b->cols;++i)
    fst(ut->data[i], &N, buf->data, &NN);
  transposeMatrix(b, ut);
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

  grid = equidistantMesh(0.0, 1.0, N);

  b = createMatrix(N-1,N-1);
  e = createMatrix(N-1,N-1);

  evalMeshInternal2(b, grid, source, local);
  h = grid->data[1]-grid->data[0];
  scaleVector(b->as_vec, pow(h, 2));
  evalMeshInternal2(e, grid, exact, local);
  axpy(b->as_vec, e->as_vec, alpha);

  lambda = generateEigenValuesP1D(N-1);

  time = WallTime();

  DiagonalizationPoisson2Dfst(b, lambda);

  printf("elapsed: %f\n", WallTime()-time);

  axpy(b->as_vec,e->as_vec,-1.0);

  printf("max error: %e\n", maxNorm(b->as_vec));

  freeMatrix(b);
  freeMatrix(e);
  freeVector(grid);
  if (lambda)
    freeVector(lambda);
  return 0;
}
