#!/bin/bash

#module add modules
module add modules/2.2-20230808
#module add modules/2.1-20230222
module add intel-oneapi-compilers
module add intel-oneapi-mpi
#module add cuda
#module add gcc
module add cuda

make clean
make -j8


#./nbodyplus.exe -f nbody.dat >stdout 2>stderr
#./nbodyplus.exe -f binary.dat >stdout 2>stderr
