#!/bin/bash

cpus="2 4 8"

for mpicpu in ${cpus}
do
    echo -e "$mpicpu"
    mpirun -np $mpicpu ex4
done
