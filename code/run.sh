#!/bin/bash
mkdir _debug
cd _debug && mpirun -np $1 ../ex4
