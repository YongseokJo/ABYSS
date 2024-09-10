export PATH=$PATH:/appl/intel/oneapi/compiler/2021.4.0/linux/bin/intel64
export PATH=/usr/local/cuda-11.5/targets/x86_64-linux/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.5/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
which icc
#module load icc/latest
make clean
make -j4 -f Makefile15
cp nbodyplus_binary.exe /data1/wispedia/nbody/purenbody
date
