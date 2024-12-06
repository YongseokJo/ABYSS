#!/bin/bash


export LIBRARY_PATH=/mnt/home/yjo10/packages/pdtoolkit-3.25.2/build/x86_64/lib:$LIBRARY_PATH
export LIBRARY_PATH=//mnt/home/yjo10/packages/tau-2.34/build/x86_64/lib:$LIBRARY_PATH


export TAU_MAKEFILE=/mnt/home/yjo10/packages/tau-2.34/build/x86_64/lib/Makefile.tau-mpi-pdt
#export TAU_MAKEFILE=/mnt/home/yjo10/packages/tau-2.34/build/x86_64/lib/Makefile.tau-mpi-pthread-pdt


#module add modules
#module add modules/2.2-20230808
#module add modules/2.1-20230222
#module add intel-oneapi-compilers
#module add intel-oneapi-mpi
module add openmpi
#module add cuda
#module add gcc
module add cuda

make clean
make -j8


#./nbodyplus.exe -f nbody.dat >stdout 2>stderr
#./nbodyplus.exe -f binary.dat >stdout 2>stderr
