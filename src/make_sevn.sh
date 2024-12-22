export LD_LIBRARY_PATH=/usr/local/cuda-11.8/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
module load icc/latest
module load mpi/latest
make clean
make sevn -j4
# cp nbodyplus_binary.exe /data1/wispedia/nbody/purenbody
date
