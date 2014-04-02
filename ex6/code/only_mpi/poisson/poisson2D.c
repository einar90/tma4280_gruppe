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
   return -30.0*pow(y,4)*x*(pow(x,5.0)-1)-30.0*pow(x,4)*y*(pow(y,5)-1);
  //return 2*x+y;
}


void safePrintMatrix(Matrix m, int proc_id)
{
  int rows = m->rows, cols = m->cols;
  char start[32], end[32];
  sprintf(start, "##   Printing from %d    ##\n", proc_id);
  sprintf(end,   "## Done printing from %d ##\n", proc_id);
  // 10 digits allocated per number. 1 digit per newline.
  int bufferLength = 10*rows*cols + 1*rows;
  int r, c, i, pos;
  char buf[bufferLength];
  pos = 0;
  for (r = 0; r < rows; ++r)
  {
    for (c = 0; c < cols; ++c)
    {
      pos += sprintf(buf + pos, "%8.4f ", m->data[c][r]);
    }
    pos += sprintf(buf + pos, "\n");
  }
  printf("%s%s%s", start, buf, end);
}


void printRawMatrix(Matrix b)
{
  int r, c;
  for (c = 0; c < b->cols; ++c)
  {
    for (r = 0; r < b->rows; ++r)
    {
      printf("%f ", b->data[c][r]);
    }
    printf("\n");
  }
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

void resetBuffer(Vector buf)
{
  int n = buf->len;
  freeVector(buf);
  buf = createVector(n);
}

void printResults(int omp, int mpi, int N, double runtime, double error)
{
  printf("omp %d\nmpi %d\nn %d\nruntime %f\nerror %e\n", omp, mpi, N, runtime, error);
}



void DiagonalizationPoisson2Dfst(Matrix b, int size, int rank,
                                 MPI_Status status, int *displacements,
                                 int *partlens)
{
  int i,j,k,r,c;
  int N=b->rows+1;
  Vector lambda;
  Matrix ut = cloneMatrix(b);
  Vector buf;
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

#pragma omp parallel for schedule(static) private(buf)
    for (i=0;i<b->cols;++i){
      buf = createVector(4*(b->rows+1));
      fst(b->data[i], &N, buf->data, &NN);
      freeVector(buf);
    }

    MPI_Send(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE, 0, 200, MPI_COMM_WORLD);
  }

  if (rank == 0)
  {
    recvParts(partlens, displacements, b, N, size, 200);
    transposeMatrix(ut, b);

    distributeCols(partlens, displacements, N, ut, copyMatrix, size, 300);
  }

  if (rank != 0)
  {
    MPI_Recv(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE,
             0, 300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for schedule(static) private(buf)
    for (i=0;i<b->cols;++i)
    {
      buf = createVector(4*(b->rows+1));
      fstinv(b->data[i], &N, buf->data, &NN);
      freeVector(buf);
    }

    MPI_Send(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE, 0, 400, MPI_COMM_WORLD);
  }

  if (rank == 0)
  {
    recvParts(partlens, displacements, ut, N, size, 400);
    for (j=0;j<b->cols;++j){
      for (i=0;i<b->rows;++i){
        ut->data[j][i] /= (lambda->data[i]+lambda->data[j]+alpha);
      }
    }

    distributeCols(partlens, displacements, N, ut, copyMatrix, size, 500);
  }

  if (rank != 0)
  {
    MPI_Recv(&ut->as_vec->data[0], ut->as_vec->len, MPI_DOUBLE,
             0, 500, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for schedule(static) private(buf)
    for (i=0;i<b->cols;++i)
    {
      buf = createVector(4*(b->rows+1));
      fst(ut->data[i], &N, buf->data, &NN);
      freeVector(buf);
    }

    MPI_Send(&ut->as_vec->data[0], ut->as_vec->len, MPI_DOUBLE, 0, 600, MPI_COMM_WORLD);
  }

  if (rank == 0)
  {
    recvParts(partlens, displacements, ut, N, size, 600);

    transposeMatrix(b, ut);

    distributeCols(partlens, displacements, N, b, copyMatrix, size, 700);
  }

  if(rank != 0) {
    MPI_Recv(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE,
             0, 700, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#pragma omp parallel for schedule(static) private(buf)
    for (i=0;i<ut->cols;++i)
    {
      buf = createVector(4*(b->rows+1));
      fstinv(b->data[i], &N, buf->data, &NN);
      freeVector(buf);
    }

    MPI_Send(&b->as_vec->data[0], b->as_vec->len, MPI_DOUBLE, 0, 800, MPI_COMM_WORLD);
  }

  if (rank == 0)
  {
    recvParts(partlens, displacements, b, N, size, 800);
  }

  freeMatrix(ut);
}

void DiagonalizationPoisson2DfstSerial(Matrix b, const Vector lambda)
{
  int i,j;
  int N=b->rows+1;
  Matrix ut = cloneMatrix(b);
  Vector buf;
  int NN=4*N;
  double time;

#pragma omp parallel for schedule(static) private(buf)
  for (i=0;i<b->cols;++i)
  {
    buf = createVector(4*(b->rows+1));
    fst(b->data[i], &N, buf->data, &NN);
    freeVector(buf);
  }

  transposeMatrix(ut, b);

#pragma omp parallel for schedule(static) private(buf)
  for (i=0;i<ut->cols;++i)
  {
    buf = createVector(4*(b->rows+1));
    fstinv(ut->data[i], &N, buf->data, &NN);
    freeVector(buf);
  }

  time = WallTime();
  for (j=0;j<b->cols;++j){
    for (i=0;i<b->rows;++i){
      ut->data[j][i] /= (lambda->data[i]+lambda->data[j]+alpha);
    }
  }

#pragma omp parallel for schedule(static) private(buf)
  for (i=0;i<b->cols;++i)
  {
    buf = createVector(4*(b->rows+1));
    fst(ut->data[i], &N, buf->data, &NN);
    freeVector(buf);
  }

  transposeMatrix(b, ut);

#pragma omp parallel for schedule(static) private(buf)
  for (i=0;i<ut->cols;++i)
  {
    buf = createVector(4*(b->rows+1));
    fstinv(b->data[i], &N, buf->data, &NN);
    freeVector(buf);
  }

  freeMatrix(ut);
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
    return 5;
  }

  local = 0;

  // Master har ansvar for dette


  // Find a split size, exclude one processor as its master and handling send/recv
  if(size == 1) {
    printf("%s\n", "Running with one prosess, openMP active");
    grid = equidistantMesh(0.0, 1.0, N);
    b = createMatrix(N-1,N-1);
    e = createMatrix(N-1,N-1);

    evalMeshInternal2(b, grid, source, local);
    h = grid->data[1]-grid->data[0];
    scaleVector(b->as_vec, pow(h, 2));
    evalMeshInternal2(e, grid, exact, local);
    axpy(b->as_vec, e->as_vec, alpha);
    Vector lambda = generateEigenValuesP1D(N-1);

    double time = WallTime();
    DiagonalizationPoisson2DfstSerial(b, lambda);
    double runtime = WallTime() - time;

    //printRawMatrix(b); // Printing solution
    axpy(b->as_vec,e->as_vec,-1.0);  // Calculating errors
    double error = maxNorm(b->as_vec);
    //printf("max error: %e\n", maxNorm(b->as_vec));

    printResults(getMaxThreads(), size, N, runtime, error);
    freeMatrix(e);
    freeVector(grid);
  }
  else
  {

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

    if(rank == 0) {
      double runtime = WallTime() - time;
      axpy(b->as_vec,e->as_vec,-1.0);  // Calculating errors
      double error = maxNorm(b->as_vec);
      printResults(getMaxThreads(), size, N, runtime, error);

      // printf("elapsed: %f\n", WallTime()-time);
      // printRawMatrix(b); // Printing solution
      // printf("max error: %e\n", maxNorm(b->as_vec));
      freeMatrix(e);
      freeVector(grid);
    }

  }
  freeMatrix(b);

  close_app();
  return 0;
}
