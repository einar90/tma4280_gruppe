#!/usr/bin/env bash
module load intelcomp/14.0.1
module load lapack/3.3.0
CXX=icpc CC=icc FC=ifort cmake ..
