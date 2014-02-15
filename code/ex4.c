#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "common.h"

void populateVector(float* array, int n) {
	for(i = 0; i < n-1; i++) {
		array[i] = 1.0 / (i*i);
	}
}

float* createVector(float* array, int n) {
	array = malloc(n * sizeof(float));
}

void openMpTest() {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < 10; ++i)
  {
    printf("OpenMP Test. Iteration %d\n", i);
  }
}


int main(int argc, char const *argv[])
{
  openMpTest();
  return 0;
}
