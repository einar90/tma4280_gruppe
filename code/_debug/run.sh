#!/bin/bash
cmake ..
make
mpirun -np $1 ../ex4
