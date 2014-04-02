#!/bin/bash
mkdir release
cd release

module load intelcomp
module load openmpi/1.4.3-intel

CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make

cp ../jobscript.sh ./jobscript.sh
chmod 755 jobscript.sh
