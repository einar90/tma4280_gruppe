#!/bin/bash
#PBS -N team_awesome_parallell_poisson
#PBS -lnodes=3:ppn=12:default
#PBS -lwalltime=00:15:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine,compact"

for N in 128 256 512 1024 2048 4096 8192 16384
do
  OMP_NUM_THREADS=6 mpirun -npernode 2  poisson2D $N
done


# Each node: 2 processors
# Each processor: 6 cores

# -> -npernode = 2
# -> omp num threads = 6
