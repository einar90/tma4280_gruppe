#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "common.h"


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
