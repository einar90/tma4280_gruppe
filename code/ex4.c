#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "common.h"

#define PI 3.14159265358979323846


void fillVector_ex4(Vector x) {
  int i;
  for(i=1; i<= x->len; ++i){
    double random = (1.0 / (double)(i*i));
    x->data[i*x->stride-1] = random;
  }
}


double sumVector(Vector vector, int start, int stop) {
  double sum = 0;
  int i;

  for(i = start; i < stop; i++) {
    sum = sum + vector->data[i];
  }
  return sum;
}

double calcError(double sum)
{
  double S = pow(PI,2.0)/6.0;
  return fabs(S - sum);
}

int intpow(int a, int b)
{
  int result = a;
  for (int i = 1; i < b; ++i)
  {
    result = result*a;
  }
  return result;
}


int main(int argc, char** argv)
{
  double starttime, endtime;
  int rank, size, i, tag;
  MPI_Status status;
  init_app(argc, argv, &rank, &size);
  tag = 100;
  Vector vector;

  int k;
  for (k = 3; k < 15; k++){
    int n = intpow(2,k);
    int partlen = n/size;
    double sum = 0;
    if (rank == 0)
    {
      starttime = WallTime();
      vector = createVector(n);
      fillVector_ex4(vector);
      for (int i = 1; i < size; ++i)
      {
        MPI_Send(&vector->data[(i-1)*partlen], partlen, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
      }
      for (int i = 1; i < size; ++i)
      {
        double recivedsum;
        MPI_Recv(&recivedsum, 1, MPI_DOUBLE, i, 42, MPI_COMM_WORLD, &status);
        sum += recivedsum;
      }
      endtime = WallTime();
      printf("%d\t", n);
      printf("%f\t", endtime-starttime);
      printf("%f\t", sum);
      printf("%f\n", calcError(sum));
    }
    else
    {
      vector = createVector(partlen);
      MPI_Recv(vector->data, vector->len, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
      double partsum = sumVector(vector, 0, partlen);
      MPI_Send(&partsum, 1, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
    }

    freeVector(vector);
  }
  close_app();
  return 0;
}

  /////////////////
  // MPI Testing //
  /////////////////
  /*
  int rank, size, i, tag;
  MPI_Status status;
  char message[20];
  init_app(argc, argv, &rank, &size);
  tag = 100;

  if (rank == 0) {
    strcpy(message, "Hello world!");
    for (i=1; i < size; ++i)
      MPI_Send(message, 13, MPI_CHAR, i, tag, MPI_COMM_WORLD);
  } else
    MPI_Recv(message, 13, MPI_CHAR, 0, tag, MPI_COMM_WORLD, &status);

  printf("process %d: %s\n", rank, message);
  */
