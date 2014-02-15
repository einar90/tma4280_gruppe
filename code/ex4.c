#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "common.h"

#define PI 3.14159265358979323846

Vector vector;

/*
   Function for testing the OpenMP lib
 */
void openMpTest() {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < 4; ++i)
  {
    printf("OpenMP Test. Iteration %d\n", i);
  }
}

void fillVector_ex4(Vector x) {
  int i;
  for(i=1; i<=x->len; ++i){
    double random = (1.0 / (double)(i*i));
    x->data[i*x->stride-1] = random;
  }
}

void setupVector(int n) {
	vector = createVector(n);
	fillVector_ex4(vector);
}

double sumVector(int start, int stop) {
  double sum = 0;
  int i = start;
  #pragma omp parallel for schedule(static) private(start) reduction(+:sum)
  for(i; i < stop; i++) {
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
  int k;
  for (k = 3; k < 15; k++){
    int n = intpow(2,k);
    setupVector(n);
    double sum = sumVector(0, n);
    printf("n = %d\n", n);
    printf("Sum: %f\n", sum);
    printf("Error S - S_n = %f\n\n", calcError(sum));
    freeVector(vector);
  }




  return 0;
}


  ////////////////////
  // OpenMP testing //
  ////////////////////
  //openMpTest();

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

  close_app();*/
