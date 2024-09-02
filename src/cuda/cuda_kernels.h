#ifndef KERNEL_H
#define KERNEL_H
#include "../defs.h"


__global__	void initialize(CUDA_REAL* result, CUDA_REAL* diff, CUDA_REAL *magnitudes, int n, int m, int* subset);
__global__ void compute_pairwise_diff_subset(const CUDA_REAL* ptcl, CUDA_REAL* diff, int n, int m, const int* subset, int start);
__global__ void compute_magnitudes_subset(const CUDA_REAL *r2, const CUDA_REAL* diff, CUDA_REAL* magnitudes, int n, int m, int* subset, bool* neighbor2, int start);
__global__ void compute_forces_subset(const CUDA_REAL* ptcl, CUDA_REAL *diff, const CUDA_REAL* magnitudes, int n, int m, const int* subset);
__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const CUDA_REAL* r2, const CUDA_REAL* magnitudes, int n, int m, const int *subset);
__global__ void reduce_forces(const CUDA_REAL *diff, CUDA_REAL *result, int n, int m);
__global__ void print_forces_subset(CUDA_REAL* result, int m);
<<<<<<< HEAD
=======

// __global__ void compute_forces(const CUDA_REAL* __restrict__ ptcl, const CUDA_REAL* __restrict__ r2, CUDA_REAL *d_diff, int m, int n, const int* __restrict__ subset, int* __restrict__ neighbor, int* num_neighbor, int start);
__global__ void compute_forces(const CUDA_REAL* __restrict__ ptcl, const CUDA_REAL* __restrict__ r2, CUDA_REAL* __restrict__ diff, int m, int n, const int* __restrict__ subset, int* __restrict__ neighbor, int* num_neighbor, int start);
__global__ void gather_neighbor(const int*  neighbor_block, const int*  num_neighbor, int* gathered_neighbor, int  m);
void force_reduction(cublasHandle_t handle, const double3 *acc, const double3 *adot, CUDA_REAL *result, int m, int n);

>>>>>>> f688615c336c02f68febfc5c6d41376999e08ae3
void reduce_forces_cublas(cublasHandle_t handle, const CUDA_REAL *diff, CUDA_REAL *result, int n, int m);

/*
__device__ void _addition(Result &result, const Result res);
__device__ void _copy(Result &result, const Result res);
*/

#endif
