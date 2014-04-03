#!/bin/bash
#PBS -N team_awesome_parallell_poisson
#PBS -lnodes=3:ppn=12:default
#PBS -lwalltime=03:00:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine,compact"

for P in 1 2 4 6 8
do
  for T in 3 4 6 8 12
  do
    OMP_NUM_THREADS=$T mpirun -npernode $P  poisson2D 512
  done
done


# Each node: 2 processors
# Each processor: 6 cores

# -> -npernode = 2
# -> omp num threads = 6
