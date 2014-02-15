#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "common.h"

/*
	Variables used in the program
*/
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
    printf("Trying to input %f at %i\n", random, i-1);
    x->data[i*x->stride-1] = random;
  }
}

void setupVector(int n) {
	vector = createVector(n);
	fillVector_ex4(vector);
}

double sumVector(Vector vec, int start, int stop) {
  double sum = 0;
  for(;start < stop; start++) {
    printf("Number is: %f, start is: %i\n", vec->data[start*vec->stride], start);
    sum = sum + vec->data[start*vec->stride];
  }
}


int main(int argc, char** argv)
{
  if (argc < 2) {
    printf("Need at least one parameter: i\n");
    return 1;
  }

  int ii = atoi(argv[1]);
  printf("Number: %i\n", ii);
  setupVector(ii);
  double crap = sumVector(vector, 0, ii);
  printf("Sum: %f\n", crap);

  ////////////////////
  // OpenMP testing //
  ////////////////////
  //openMpTest();

  /////////////////
  // MPI Testing //
  /////////////////
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

  close_app();

  return 0;
}
