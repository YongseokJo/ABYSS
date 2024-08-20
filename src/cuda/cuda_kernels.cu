#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cassert>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "cuda_defs.h"
#include "../defs.h"
#include "cuda_kernels.h"

#ifdef THRUST
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/find.h>
#endif


// CUDA kernel to compute the forces for a subset of particles
__global__ void print_forces_subset(CUDA_REAL* result, int m) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m) {
		printf("acc_x: (%d) = %e\n", idx, result[_six*idx]);
		printf("acc_y: (%d) = %e\n", idx, result[_six*idx+1]);
		printf("acc_z: (%d) = %e\n", idx, result[_six*idx+2]);
		printf("adot_x: (%d) = %e\n", idx, result[_six*idx+3]);
		printf("adot_y: (%d) = %e\n", idx, result[_six*idx+4]);
		printf("adot_z: (%d) = %e\n", idx, result[_six*idx+5]);
				/*
				atomicAdd(&result[i+1], scale * diff[six_idx + 1]);
				atomicAdd(&result[i+2], scale * diff[six_idx + 2]);

				atomicAdd(&result[i+3], scale * (diff[six_idx + 3] - magnitudes[idx+1]*diff[six_idx    ]/magnitudes[idx]));
				atomicAdd(&result[i+4], scale * (diff[six_idx + 4] - magnitudes[idx+1]*diff[six_idx + 1]/magnitudes[idx]));
				atomicAdd(&result[i+5], scale * (diff[six_idx + 5] - magnitudes[idx+1]*diff[six_idx + 2]/magnitudes[idx]));
				*/
	}
}

// NTHREAED = 64;
// NJBlock = 28
// NNB_per_block = 256;
__global__ void compute_forces(const CUDA_REAL* __restrict__ ptcl, const CUDA_REAL* __restrict__ r2, CUDA_REAL* __restrict__ diff, int m, int n, const int* __restrict__ subset, int* __restrict__ neighbor, int* num_neighbor, int start){
	// define i and j. in this code, grid is 2D and block is 1D
    int i = threadIdx.x + blockIdx.x * blockDim.x; // Unique thread index across all blocks
	int tid = threadIdx.x;
	// int BatchSize = blockDim.x;
	// int j_begin = blockIdx.x * (n / (BatchSize * gridDim.x));
	int idx_save_size = gridDim.y * m;

	int j_begin = blockIdx.y * n / gridDim.y;
	int j_end = (blockIdx.y + 1) * n / gridDim.y;
	if (blockIdx.y == gridDim.y - 1) j_end = n;  // Ensure the last block covers all remaining elements
	
	while (i < m + BatchSize - 1){
		int i_ptcl = subset[i + start];
		if (i >= m) i_ptcl = subset[m - 1 + start]; //assign dummy values for the last block	
		CUDA_REAL pi_x = ptcl[i_ptcl];
		CUDA_REAL pi_y = ptcl[i_ptcl + n];
		CUDA_REAL pi_z = ptcl[i_ptcl + 2 * n];
		CUDA_REAL pi_vx = ptcl[i_ptcl + 3 * n];
		CUDA_REAL pi_vy = ptcl[i_ptcl + 4 * n];
		CUDA_REAL pi_vz = ptcl[i_ptcl + 5 * n];
		int NumNeighbor = 0;
		int idx_save = i * gridDim.y + blockIdx.y;
		double save_acc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		int* BlockNeighbor = &neighbor[NNB_per_block*idx_save]; // Pointer to the neighbor list of the current block

		for (int j=j_begin; j < j_end; j+=BatchSize){
			int current_batch_size = min(BatchSize, n - j);

			// assing shared particles for BatchSize particles to each block
			// __shared__ CUDA_REAL sh_ptcl[BatchSize*7];
			__shared__ double sh_pos_x[BatchSize];
			__shared__ double sh_pos_y[BatchSize];
			__shared__ double sh_pos_z[BatchSize];
			__shared__ double sh_vel_x[BatchSize];
			__shared__ double sh_vel_y[BatchSize];
			__shared__ double sh_vel_z[BatchSize];
			__shared__ double sh_mass[BatchSize];
			__syncthreads();
			if (tid < current_batch_size) {
				sh_pos_x[tid] = ptcl[j + tid];
				sh_pos_y[tid] = ptcl[j + tid + n];
				sh_pos_z[tid] = ptcl[j + tid + 2 * n];
				sh_vel_x[tid] = ptcl[j + tid + 3 * n];
				sh_vel_y[tid] = ptcl[j + tid + 4 * n];
				sh_vel_z[tid] = ptcl[j + tid + 5 * n];
				sh_mass[tid]  = ptcl[j + tid + 6 * n];
			}
			__syncthreads();

            #pragma unroll 4
			for (int jj=0; jj<current_batch_size; jj++){

				if (i<m){
					// Calculate forces
					CUDA_REAL dx = sh_pos_x[jj] - pi_x;
					CUDA_REAL dy = sh_pos_y[jj] - pi_y;
					CUDA_REAL dz = sh_pos_z[jj] - pi_z;
					CUDA_REAL magnitude0 = dx*dx + dy*dy + dz*dz;
					// int idx = i * n + (j + jj);
					// neighbor[idx] = isNeighbor;

					if (magnitude0 > r2[i_ptcl]) {
						// Calculate velocity differences
						CUDA_REAL dvx = sh_vel_x[jj] - pi_vx;
						CUDA_REAL dvy = sh_vel_y[jj] - pi_vy;
						CUDA_REAL dvz = sh_vel_z[jj] - pi_vz;
						CUDA_REAL inv_sqrt_m0 = rsqrt(magnitude0);
						CUDA_REAL inv_m0 = inv_sqrt_m0 * inv_sqrt_m0; // or 1 / magnitude0
						CUDA_REAL scale = sh_mass[jj] * inv_sqrt_m0 * inv_m0;
						// Calculate adot_temp
						CUDA_REAL common_factor = 3.0 * (dx*dvx + dy*dvy + dz*dvz) * inv_m0;
						
						save_acc[0] += scale * dx;
						save_acc[1] += scale * dy;
						save_acc[2] += scale * dz;
						save_acc[3] += scale * (dvx - common_factor * dx);
						save_acc[4] += scale * (dvy - common_factor * dy);
						save_acc[5] += scale * (dvz - common_factor * dz);

					}
					else if (i_ptcl != j+jj) {
						BlockNeighbor[NumNeighbor++] = j + jj;
						assert (NumNeighbor < NNB_per_block);
					}
				}
			} //end of jj loop
		}//end of j loop
		if (i < m){
			diff[idx_save] = save_acc[0];
			diff[idx_save + idx_save_size] = save_acc[1];
			diff[idx_save + 2 * idx_save_size] = save_acc[2];
			diff[idx_save + 3 * idx_save_size] = save_acc[3];
			diff[idx_save + 4 * idx_save_size] = save_acc[4];
			diff[idx_save + 5 * idx_save_size] = save_acc[5];
			num_neighbor[idx_save] = NumNeighbor;
		}
		i += gridDim.x * blockDim.x;
	} //end of i loop
}

// new kernel to compute the forces
__global__ void compute_forces_old(const CUDA_REAL* __restrict__ ptcl, const CUDA_REAL* __restrict__ r2, CUDA_REAL* __restrict__ diff, int m, int n, const int* __restrict__ subset, bool* __restrict__ neighbor, int start){
	int j = blockIdx.x * blockDim.x + threadIdx.x; // Index for particles (j)
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Index for subset (i)

	if (i < m && j < n) {
		int i_ptcl = subset[i + start];
		int idx = i * n + j;
		
        // Load particle data into registers for better memory access efficiency
        CUDA_REAL pj_x = ptcl[j];
        CUDA_REAL pj_y = ptcl[j + n];
        CUDA_REAL pj_z = ptcl[j + 2 * n];
        CUDA_REAL pj_vx = ptcl[j + 3 * n];
        CUDA_REAL pj_vy = ptcl[j + 4 * n];
        CUDA_REAL pj_vz = ptcl[j + 5 * n];
        CUDA_REAL pj_mass = ptcl[j + 6 * n];

        CUDA_REAL pi_x = ptcl[i_ptcl];
        CUDA_REAL pi_y = ptcl[i_ptcl + n];
        CUDA_REAL pi_z = ptcl[i_ptcl + 2 * n];
        CUDA_REAL pi_vx = ptcl[i_ptcl + 3 * n];
        CUDA_REAL pi_vy = ptcl[i_ptcl + 4 * n];
        CUDA_REAL pi_vz = ptcl[i_ptcl + 5 * n];

        // Calculate position differences
        CUDA_REAL dx = pj_x - pi_x;
        CUDA_REAL dy = pj_y - pi_y;
        CUDA_REAL dz = pj_z - pi_z;

        // CUDA_REAL magnitude0 = dx*dx + dy*dy + dz*dz;
		CUDA_REAL magnitude0 = dx*dx + dy*dy + dz*dz;

        // Initialize result variables
        double ax = 0.0, ay = 0.0, az = 0.0;
        double adotx = 0.0, adoty = 0.0, adotz = 0.0;

        bool isNeighbor = magnitude0 <= r2[i_ptcl];
		neighbor[idx] = isNeighbor;

        if (!isNeighbor) {
            // Calculate velocity differences
            CUDA_REAL dvx = pj_vx - pi_vx;
            CUDA_REAL dvy = pj_vy - pi_vy;
            CUDA_REAL dvz = pj_vz - pi_vz;

            // CUDA_REAL magnitude1 = dx*dvx + dy*dvy + dz*dvz;
			CUDA_REAL magnitude1 = dx*dvx + dy*dvy + dz*dvz;
			CUDA_REAL inv_sqrt_m0 = rsqrt(magnitude0);
			CUDA_REAL inv_m0 = inv_sqrt_m0 * inv_sqrt_m0; // or 1 / magnitude0

            CUDA_REAL scale = ptcl[j + 6 * n] * inv_sqrt_m0 * inv_m0;
            ax = scale * dx;
            ay = scale * dy;
            az = scale * dz;

            // Calculate adot_temp
            CUDA_REAL common_factor = 3.0 * magnitude1 * inv_m0;
            adotx = scale * (dvx - common_factor * dx);
            adoty = scale * (dvy - common_factor * dy);
            adotz = scale * (dvz - common_factor * dz);
	        }

        // Store results in the diff array
        diff[idx] = ax;
        diff[idx + n * m] = ay;
        diff[idx + 2 * n * m] = az;
        diff[idx + 3 * n * m] = adotx;
        diff[idx + 4 * n * m] = adoty;
        diff[idx + 5 * n * m] = adotz;
		//printf("a0 (%d,%d) = %e, %e %e\n \t a1 = %e, %e %e\n", i, j, diff[six_idx], diff[six_idx+1], diff[six_idx+2], diff[six_idx+3], diff[six_idx+4], diff[six_idx+5]);

    }
}

void force_reduction(cublasHandle_t handle, const double3 *acc, const double3 *adot, CUDA_REAL *result, int m, int n) {
    CUDA_REAL *d_matrix;
    cudaMalloc(&d_matrix, m * n * sizeof(CUDA_REAL));
	
    // Create a vector of ones for the summation
    double *ones;
    cudaMalloc(&ones, n * sizeof(double));
    double *h_ones = new double[n];
    for (int i = 0; i < n; ++i) {
        h_ones[i] = 1.0;
    }
    cudaMemcpy(ones, h_ones, n * sizeof(double), cudaMemcpyHostToDevice);

    // Initialize result array to zero
    cudaMemset(result, 0, m * 6 * sizeof(double));

    const double alpha = 1.0;
    const double beta = 0.0;

    // Sum over the second axis (n) for each of the 6 elements (acc + adot)
    for (int i = 0; i < 3; ++i) {
        // Copy acc components to the matrix
        cublasDcopy(handle, m * n, reinterpret_cast<const double*>(acc) + i, 3, d_matrix, 1);
        cublasDgemv(
            handle,
            CUBLAS_OP_T,  // Transpose
            n,            // Number of rows of the matrix A
            m,            // Number of columns of the matrix A
            &alpha,       // Scalar alpha
            d_matrix,     // Pointer to the first element of the i-th sub-matrix
            n,            // Leading dimension of the sub-matrix
            ones,         // Pointer to the vector x
            1,            // Increment between elements of x
            &beta,        // Scalar beta
            result + i,   // Pointer to the first element of the result vector
            6             // Increment between elements of the result vector
        );

        // Copy adot components to the matrix
        cublasDcopy(handle, m * n, reinterpret_cast<const double*>(adot) + i, 3, d_matrix, 1);
        cublasDgemv(
            handle,
            CUBLAS_OP_T,  // Transpose
            n,            // Number of rows of the matrix A
            m,            // Number of columns of the matrix A
            &alpha,       // Scalar alpha
            d_matrix,     // Pointer to the first element of the i-th sub-matrix
            n,            // Leading dimension of the sub-matrix
            ones,         // Pointer to the vector x
            1,            // Increment between elements of x
            &beta,        // Scalar beta
            result + 3 + i, // Pointer to the first element of the result vector
            6             // Increment between elements of the result vector
        );
    }
    // Cleanup
    delete[] h_ones;
    cudaFree(ones);
    cudaFree(d_matrix);
}



__global__	void initialize(CUDA_REAL* result, CUDA_REAL* diff, int n, int m, int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < m*n) {
		diff[_six*idx    ] = 0.;
		diff[_six*idx + 1] = 0.;
		diff[_six*idx + 2] = 0.;
		diff[_six*idx + 3] = 0.;
		diff[_six*idx + 4] = 0.;
		diff[_six*idx + 5] = 0.;
	}

	#ifdef old
	if (idx < m * n) {
		int i = idx / n;
		int j = idx % n;

		if (j == 0) {
			result[_six*i] = 0.;
			result[_six*i + 1] = 0.;
			result[_six*i + 2] = 0.;
			result[_six*i + 3] = 0.;
			result[_six*i + 4] = 0.;
			result[_six*i + 5] = 0.;
			// num_neighbor[i] = 0;
			/*
			for (j=0; j<NumNeighborMax; j++)
				neighbor[NumNeighborMax*i+j] = 0;
				*/
		}
	}
	#endif
}

// CUDA kernel to compute pairwise differences for a subset of particles
__global__ void compute_pairwise_diff_subset(const CUDA_REAL* ptcl, CUDA_REAL* diff, int n, int m, const int* subset, int start) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		int i = subset[idx / n + start];
		int j = idx % n;
		idx *= _six;
		i *= _seven;
		j *= _seven;

		diff[idx]   = ptcl[j]   - ptcl[i];
		diff[idx+1] = ptcl[j+1] - ptcl[i+1];
		diff[idx+2] = ptcl[j+2] - ptcl[i+2];
		diff[idx+3] = ptcl[j+3] - ptcl[i+3];
		diff[idx+4] = ptcl[j+4] - ptcl[i+4];
		diff[idx+5] = ptcl[j+5] - ptcl[i+5];

		//printf("(%d,%d) = %e, %e, %e\n", i/_seven, j/_seven,  ptcl[i], ptcl[j], diff[idx]);
	}
}

// n: NNB, m: NumTarget
__global__ void compute_magnitudes_subset(const CUDA_REAL *r2, const CUDA_REAL* diff, CUDA_REAL* magnitudes, int n, int m, int* subset, bool* neighbor2, int start) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < n * m) {
		int i = subset[idx / n + start];
		int j = idx % n;
		int six_idx = _six*idx;
		int two_idx = _two*idx;


		magnitudes[two_idx]   += diff[(six_idx)]    *diff[(six_idx)];
		magnitudes[two_idx]   += diff[(six_idx) + 1]*diff[(six_idx) + 1];
		magnitudes[two_idx]   += diff[(six_idx) + 2]*diff[(six_idx) + 2];
		magnitudes[two_idx+1] += diff[(six_idx)]    *diff[(six_idx) + 3];
		magnitudes[two_idx+1] += diff[(six_idx) + 1]*diff[(six_idx) + 4];
		magnitudes[two_idx+1] += diff[(six_idx) + 2]*diff[(six_idx) + 5];

		//printf("(%d,%d) = %e, %e\n", i, j,  magnitudes[two_idx], r2[i]);

		if (magnitudes[two_idx] <= r2[i]) {
			//printf("(%d, %d, %d): %e, %e\n", idx/n, i, j, magnitudes[two_idx], r2[i]);
			magnitudes[two_idx]   = -magnitudes[two_idx];
			neighbor2[idx] = true;
		}
		else {
			neighbor2[idx] = false;
		}
	}
}


// CUDA kernel to compute the forces for a subset of particles
__global__ void compute_forces_subset(const CUDA_REAL* ptcl, CUDA_REAL *diff, const CUDA_REAL* magnitudes, int n, int m, const int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		//int i = subset[idx / n];
		int i = idx / n;
		int j = idx % n;
		int six_idx = idx*_six;
		idx *= _two;
		CUDA_REAL acc[Dim], adot[Dim];

		if (magnitudes[idx] <= 0.) {
			acc[0]  = 0.;
			acc[1]  = 0.;
			acc[2]  = 0.;
			adot[0] = 0.;
			adot[1] = 0.;
			adot[2] = 0.;
		}
		else {
			CUDA_REAL scale = ptcl[_seven*j+6] / (magnitudes[idx] * sqrtf(magnitudes[idx]));
			acc[0]  = scale * diff[six_idx];
			acc[1]  = scale * diff[six_idx + 1];
			acc[2]  = scale * diff[six_idx + 2];

			adot[0] = scale * (diff[six_idx + 3] - 3*magnitudes[idx+1]*diff[six_idx    ]/magnitudes[idx]);
			adot[1] = scale * (diff[six_idx + 4] - 3*magnitudes[idx+1]*diff[six_idx + 1]/magnitudes[idx]);
			adot[2] = scale * (diff[six_idx + 5] - 3*magnitudes[idx+1]*diff[six_idx + 2]/magnitudes[idx]);
		}

		diff[six_idx]   = acc[0];
		diff[six_idx+1] = acc[1];
		diff[six_idx+2] = acc[2];
		diff[six_idx+3] = adot[0];
		diff[six_idx+4] = adot[1];
		diff[six_idx+5] = adot[2];

		//printf("compute_forces: (%d, %d) = %e\n", i, j,  diff[six_idx]);
	}
}



/*
__device__ CUDA_REAL warpReduce(CUDA_REAL val) {
	val += __shfl_down_sync(0xffffffff, val, 16);
	val += __shfl_down_sync(0xffffffff, val, 8);
	val += __shfl_down_sync(0xffffffff, val, 4);
	val += __shfl_down_sync(0xffffffff, val, 2);
	val += __shfl_down_sync(0xffffffff, val, 1);
	return val;
}
*/

__inline__ __device__ CUDA_REAL warpReduce(CUDA_REAL val)
{
	for (int offset = warpSize/2; offset > 0; offset /= 2) 
		val += __shfl_down_sync(0xffffffff, val, offset);
	return val;
}


#define NEW_FORCE
#ifdef NEW_FORCE
__global__ void reduce_forces(const CUDA_REAL *diff, CUDA_REAL *result, int n, int m) {

	__shared__ CUDA_REAL warpSum[64]; // Assumes max 32 warps per block
	__shared__ CUDA_REAL res[_six]; //  this is for storing the results
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;
	int bdim = blockDim.x;
	CUDA_REAL sum;
	int six_idx;
	int k,l, i = blockIdx.x, j;
	int a = (n+bdim-1)/(bdim);

	if (threadIdx.x < _six) 
		res[threadIdx.x] = 0.;
	__syncthreads();


	for (l=0;l<a*bdim;l+=bdim) {
		j = threadIdx.x + l;
		//printf("(%d,%d,%d)\n", i, j, l);
		six_idx = _six*(i*n+j);
		warpSum[wid] = 0.;
		__syncthreads();
#pragma unroll 
		for (k=0;k<_six;k++) { // ax ay az adotx adoty adotz

			sum = (i < m && j < n) ? diff[six_idx+k] : 0;

			/*
			if (k == 0)
				if (i < m && j < n)
					printf("(%d,%d) = %e, %e\n", i, j, sum, diff[six_idx+k]);
					*/

			// Warp reduce
			sum = warpReduce(sum);
			/*
			if (k == 0)
				if (i < m && j < n)
					printf("first reduction (%d,%d) = %e\n", i, j, sum);
					*/

			// Block reduce
			if (lane == 0) warpSum[wid] = sum;
			__syncthreads();

			if (wid == 0)
			{
				sum = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane] : 0;
				/*
				if (k == 0)
					if (i < m && j < n)
						printf("before second reduction (%d,%d) = %e\n", i, j, sum);
						*/

				sum = warpReduce(sum);

				/*
				if (k == 0)
					if (i < m && j < n)
						printf("second reduction (%d,%d) = %e\n", i, j, sum);
						*/

				if (lane == 0 && i < m) {
					res[k] += sum;
					//printf("%d = (%e,%e)\n", i, sum, res[k]);
				}
			}
			__syncthreads();
		} // reduce across threads
	}
	if (wid == 0 && lane == 0 && i < m) {
#pragma unroll
		for (k=0; k<_six;k++) {
			//printf("%d = (%e)\n", threadIdx.x, res[k]);
			result[_six*i+k] = res[k];
		}
	}
	__syncthreads();
}

#else
__global__ void reduce_forces(const CUDA_REAL *diff, CUDA_REAL *result, int n, int m) {
	int idx = blockIdx.x * n + threadIdx.x;
	__shared__ CUDA_REAL warpSum[64]; // Assumes max 32 warps per block
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;
	CUDA_REAL sum;
	int i = blockIdx.x;
	int j = threadIdx.x;
	int six_idx = _six*(i*n+j);
	int k;

	//printf("old version\n");
#pragma unroll 
	for (k=0;k<_six;k++) {
		sum = (i < m && j < n) ? diff[six_idx+k] : 0;
		/*
		if (k == 0)
			if (i < m && j < n)
				printf("(%d,%d) = %e\n", blockIdx.x,threadIdx.x, diff[six_idx+k]);
				*/

		// Warp reduce
		sum = warpReduce(sum);

		// Block reduce
		if (lane == 0) warpSum[wid] = sum;
		__syncthreads();

		if (wid == 0)
		{
			sum = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane] : 0;
			sum = warpReduce(sum);
			if (lane == 0) result[_six*i+k] = sum;
		}
	}
}
#endif




#define MAX_SIZE 9 // maximum size of int array  blockDim.x*MaxSize = total size of int array 
#define NEW_V2 
#ifdef NEW_V2 // this works fine as :)
__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const CUDA_REAL* r2, const CUDA_REAL* magnitudes, int n, int m, const int *subset) {

	//int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	int bid = blockIdx.x;

		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bid,n,blockDim.x);
	if (tid < n && bid < m) {
		extern __shared__ int sdata[];
		__shared__ int offset;

		sdata[tid+1]=0;
		if (tid==0)  {
			sdata[tid]=0;
			offset = 0;
		}
		__syncthreads();

		int gdim = min(m,gridDim.x);
		int bdim = min(n,blockDim.x);
		int list[MAX_SIZE];
		int n_num;
		int a = (n+bdim*MAX_SIZE-1)/(bdim*MAX_SIZE);
		int start, end;
		int i = subset[bid]; //target particle id
		int j, k;
		int idx = 0;

		//printf("assign_neighbor: %d\n",l);

		// background particles I
		for (k=0; k<a; k++) {

			start = tid+(k*MAX_SIZE*bdim);
			end   = min((k+1)*MAX_SIZE*bdim, n);
			sdata[tid+1] = 0;
			n_num = 0;

			//printf("tid=%d: (start, end) =(%d, %d)\n", tid, start, end);

			// background particles II
			for (j=start; j<end; j+=bdim) {
				if (i != j) {
					//printf("(l,j)=(%d,%d)\n",l,j);
					idx = _two*(n*bid+j);
					//printf("(%d, %d,%d) = %d, %e, %e\n", l, i, j, num_neighbor[l], magnitudes[idx], r2[l]);
					if (magnitudes[idx] < 0) {
						list[n_num] = j;
						n_num++;
						//printf("(%d,%d,%d) = %d, %e, %e\n", l, i, j, n_num, magnitudes[idx], r2[l]);
					}
				}
			} // endfor bk ptcl II
			sdata[tid+1] = n_num;

			__syncthreads();

			if (tid == 0) {

				for (j=2; j<=bdim; j++)
					sdata[j] += sdata[j-1];

				if ((offset+sdata[bdim]) > NumNeighborMax) {
					printf("blockid=%d, Too many neighbors (%d, %d)\n", bid, offset, sdata[bdim]);
					assert(offset+sdata[bdim] < NumNeighborMax);
				}

				/*
					 if (l == 11) {
					 printf("\n(%d, %d) = sdata[bdim]=%d\n", l, i, sdata[bdim]);
				//for (j=0; j<=bdim; j++) 
				//printf("%d, ",sdata[j]);
				//printf("\n");
				}
				 */
			}
			__syncthreads();

			/*
				 if (l==0)
				 printf("(%d,%d), (num, sdata) =%d, %d\n", l, tid, n_num, sdata[tid]);
			 */

			for (j=0;j<n_num;j++) {
				neighbor[NumNeighborMax*bid+offset+sdata[tid]+j] = list[j];
				//printf("(%d,%d), j=%d\n", l, tid, list[j]);
			}
			__syncthreads();

			if (tid == 0) {
				//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
				offset += sdata[bdim];
				//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
			}
			__syncthreads();
		} //endfor bk ptcl I

		if (tid == 0) {
			num_neighbor[bid] = offset; // bid shoud be modified
			offset = 0;
			sdata[0] = 0;
		}
	} // m*n stuff
}

#elif defined(NEW_V1) // works well

__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const CUDA_REAL* r2, const CUDA_REAL* magnitudes, int n, int m, const int *subset) {

	//int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	int bid = blockIdx.x;

		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bid,n,blockDim.x);
	if (tid < n && bid < m) {
		extern __shared__ int sdata[];
		__shared__ int offset;

		sdata[tid+1]=0;
		if (tid==0)  {
			sdata[tid]=0;
			offset = 0;
		}
		__syncthreads();

		int gdim = min(m,gridDim.x);
		int bdim = min(n,blockDim.x);
		int list[MAX_SIZE];
		int n_num;
		int a = (n+bdim*MAX_SIZE-1)/(bdim*MAX_SIZE);
		int start, end;
		int i; //target particle id
		int j, k, l;
		int idx = 0;


		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bdim,n,blockDim.x);
		// target particles
		for (l=bid; l<m; l+=gdim) {
			i = subset[l];
			//printf("assign_neighbor: %d\n",l);

			// background particles I
			for (k=0; k<a; k++) {

				start = tid+(k*MAX_SIZE*bdim);
				end   = min((k+1)*MAX_SIZE*bdim, n);
				sdata[tid+1] = 0;
				n_num = 0;

				//printf("tid=%d: (start, end) =(%d, %d)\n", tid, start, end);

				// background particles II
				for (j=start; j<end; j+=bdim) {
					if (i != j) {
						//printf("(l,j)=(%d,%d)\n",l,j);
						idx = _two*(n*l+j);
						//printf("(%d, %d,%d) = %d, %e, %e\n", l, i, j, num_neighbor[l], magnitudes[idx], r2[l]);
						if (magnitudes[idx] < 0) {
							list[n_num] = j;
							n_num++;
							//printf("(%d,%d,%d) = %d, %e, %e\n", l, i, j, n_num, magnitudes[idx], r2[l]);
						}
					}
				} // endfor bk ptcl II
				sdata[tid+1] = n_num;

				__syncthreads();

				if (tid == 0) {

					for (j=2; j<=bdim; j++)
						sdata[j] += sdata[j-1];

					if ((offset+sdata[bdim]) > NumNeighborMax) {
						printf("blockid=%d, Too many neighbors (%d, %d)\n", bid, offset, sdata[bdim]);
						assert(offset+sdata[bdim] < NumNeighborMax);
					}

					/*
					if (l == 11) {
						printf("\n(%d, %d) = sdata[bdim]=%d\n", l, i, sdata[bdim]);
						//for (j=0; j<=bdim; j++) 
							//printf("%d, ",sdata[j]);
						//printf("\n");
					}
					*/
				}
				__syncthreads();

				/*
					 if (l==0)
					 printf("(%d,%d), (num, sdata) =%d, %d\n", l, tid, n_num, sdata[tid]);
				 */

				for (j=0;j<n_num;j++) {
					neighbor[NumNeighborMax*l+offset+sdata[tid]+j] = list[j];
					//printf("(%d,%d), j=%d\n", l, tid, list[j]);
				}
				__syncthreads();

				if (tid == 0) {
					//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
					offset += sdata[bdim];
					//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
				}
				__syncthreads();
			} //endfor bk ptcl I
			if (tid == 0) {
				num_neighbor[l] = offset; // bid shoud be modified
				offset = 0;
				sdata[0] = 0;
			}
			__syncthreads();
		} // endfor target paticles
	} // m*n stuff
}





#else

__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const REAL* r2, const REAL* magnitudes, int n, int m, const int *subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m) {
		int i = subset[idx];
		int k = 0;

		for (int j = 0; j < n; j++) {
			if (i != j) {
				k = _two*(n*idx+j);
				if (magnitudes[k] < 0) {
					//printf("(%d, %d,%d) = %d, %e, %e\n", idx, i, j, num_neighbor[idx], magnitudes[k], r2[i]);
					neighbor[NumNeighborMax*idx+num_neighbor[idx]] = j;
					num_neighbor[idx]++;
					if (num_neighbor[idx] > 100)  {
						//printf("Error: (%d, %d,%d) = %d, %e, %e\n", idx, i, j, num_neighbor[idx], magnitudes[k], r2[i]);
						assert(num_neighbor[idx] < 100);
						return;
					}
				}
			}
		}
	}
}

#endif

void reduce_forces_cublas(cublasHandle_t handle, const CUDA_REAL *diff, CUDA_REAL *result, int n, int m) {

	// CUDA_REAL *d_matrix;
    // cudaMalloc(&d_matrix, m * n * sizeof(CUDA_REAL));

    // Create a vector of ones for the summation
    double *ones;
    cudaMalloc(&ones, n * sizeof(double));
    double *h_ones = new double[n];
    for (int i = 0; i < n; ++i) {
        h_ones[i] = 1.0;
    }
    cudaMemcpy(ones, h_ones, n * sizeof(double), cudaMemcpyHostToDevice);
    // Initialize result array to zero
    // cudaMemset(result, 0, m * 6 * sizeof(double));

    const double alpha = 1.0;
    const double beta = 0.0;

    // Sum over the second axis (n) for each of the 6 elements
    for (int i = 0; i < _six; ++i) {

		// cublasDcopy(handle, m * n, diff + i, _six, d_matrix, 1);
		const CUDA_REAL *component_diff = diff + i * m * n;

        cublasDgemv(
            handle,
            CUBLAS_OP_T,  // Transpose
            n,            // Number of rows of the matrix A
            m,            // Number of columns of the matrix A
            &alpha,       // Scalar alpha
            component_diff, // Pointer to the first element of the i-th sub-matrix
            n,     // Leading dimension of the sub-matrix
            ones,         // Pointer to the vector x
            1,            // Increment between elements of x
            &beta,        // Scalar beta
            result + i, // Pointer to the first element of the result vector
            _six             // Increment between elements of the result vector
        );
    }
    // Cleanup
    delete[] h_ones;
    cudaFree(ones);
	// cudaFree(d_matrix);
}

#ifdef THRUST

struct less_than_zero
{
    __host__ __device__ bool operator()(const float x) const
    {
        return x < 0;
    }
};


void reduce_forces_thrust(const CUDA_REAL *diff, CUDA_REAL *result, int n, int m) {
    // Wrap raw pointers with Thrust device pointers
    thrust::device_ptr<const CUDA_REAL> d_diff(diff);
    thrust::device_ptr<CUDA_REAL> d_result(result);

    // Initialize result array to zero
    thrust::fill(d_result, d_result + m * 6, 0);

    // Sum over the second axis (n) for each of the 6 elements
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < m; ++j) {
            // Calculate the start and end pointers for the current sub-matrix
            thrust::device_ptr<const CUDA_REAL> start = d_diff + i + j * n * 6;
            thrust::device_ptr<const CUDA_REAL> end = start + n * 6;

            // Create a thrust device vector from start to end
            thrust::device_vector<CUDA_REAL> sub_matrix(start, end);

            // Reduce the sub-matrix and store the result
            d_result[i + j * 6] = thrust::reduce(sub_matrix.begin(), sub_matrix.end());
        }
    }
}


void reduce_neighbors(cublasHandle_t handle, int *neighbor, int* num_neighbor, CUDA_REAL *magnitudes, int n, int m, int* subset) {

	CUDA_REAL *d_matrix;
    cudaMalloc(&d_matrix, m * n * sizeof(CUDA_REAL));
    cublasDcopy(handle, m * n, magnitudes, _two, d_matrix, 1);


	for (int row = 0; row < m; ++row){
		CUDA_REAL val = 1.0;
        cudaMemcpy(d_matrix + row * n + subset[row], &val, sizeof(CUDA_REAL), cudaMemcpyHostToDevice);
	}

    // Wrap raw device pointers with thrust device pointers
    thrust::device_ptr<const CUDA_REAL> d_ptr(d_matrix);
    thrust::device_ptr<int> d_neighbor(neighbor);
    thrust::device_ptr<int> d_num_neighbor(num_neighbor);

    // Process each row
    for (int row = 0; row < m; ++row) {
        auto row_start = d_ptr + row * n;
        auto row_end = row_start + n;

        thrust::counting_iterator<int> index_sequence(0);

        // Use thrust::copy_if to select indices where elements are less than zero
        auto end = thrust::copy_if(index_sequence, index_sequence + n, row_start, d_neighbor + row * NumNeighborMax, less_than_zero());

        // Calculate the number of negative elements in the current row
        int num_neg_elements = thrust::distance(d_neighbor + row * NumNeighborMax, end);

        if (num_neg_elements > NumNeighborMax) {
            cudaFree(d_matrix);
            throw std::runtime_error("Number of negative elements exceeds NumNeighborMax");
        }

        d_num_neighbor[row] = num_neg_elements;
    }

    cudaFree(d_matrix);
}


__global__ void compute_forces_test2(const CUDA_REAL* __restrict__ ptcl, const CUDA_REAL* __restrict__ r2, CUDA_REAL* __restrict__ diff, int m, int n, const int* __restrict__ subset, bool* __restrict__ neighbor, int start){
    extern __shared__ CUDA_REAL shared_ptcl[]; // Shared memory for particle data

    int tid = threadIdx.x; // Thread ID within the block
    int i = blockIdx.y * blockDim.y + threadIdx.y; // Index for subset (i)
    int tile_size = blockDim.x; // Number of threads per block
    int tiles = (n + tile_size - 1) / tile_size; // Number of tiles required to cover all particles

    // Initialize shared memory for storing particle data
    CUDA_REAL* sh_x = shared_ptcl;
    CUDA_REAL* sh_y = sh_x + tile_size;
    CUDA_REAL* sh_z = sh_y + tile_size;
    CUDA_REAL* sh_vx = sh_z + tile_size;
    CUDA_REAL* sh_vy = sh_vx + tile_size;
    CUDA_REAL* sh_vz = sh_vy + tile_size;
    CUDA_REAL* sh_mass = sh_vz + tile_size;
	
	if (i < m) {
        int i_ptcl = subset[i + start];
        CUDA_REAL pi_x = ptcl[i_ptcl];
        CUDA_REAL pi_y = ptcl[i_ptcl + n];
        CUDA_REAL pi_z = ptcl[i_ptcl + 2 * n];
        CUDA_REAL pi_vx = ptcl[i_ptcl + 3 * n];
        CUDA_REAL pi_vy = ptcl[i_ptcl + 4 * n];
        CUDA_REAL pi_vz = ptcl[i_ptcl + 5 * n];

        for (int t = 0; t < tiles; ++t) {
            int j = t * tile_size + tid; // Global index for particle j
            if (j < n) {
                // Load particle data into shared memory
                sh_x[tid] = ptcl[j];
                sh_y[tid] = ptcl[j + n];
                sh_z[tid] = ptcl[j + 2 * n];
                sh_vx[tid] = ptcl[j + 3 * n];
                sh_vy[tid] = ptcl[j + 4 * n];
                sh_vz[tid] = ptcl[j + 5 * n];
                sh_mass[tid] = ptcl[j + 6 * n];
            }
            __syncthreads();

            if (j < n) {
                int idx = i * n + j;

                // Calculate position differences
                CUDA_REAL dx = sh_x[tid] - pi_x;
                CUDA_REAL dy = sh_y[tid] - pi_y;
                CUDA_REAL dz = sh_z[tid] - pi_z;

                CUDA_REAL magnitude0 = dx * dx + dy * dy + dz * dz;

                // Initialize result variables
                double ax = 0.0, ay = 0.0, az = 0.0;
                double adotx = 0.0, adoty = 0.0, adotz = 0.0;

                bool isNeighbor = magnitude0 <= r2[i_ptcl];
                neighbor[idx] = isNeighbor;

                if (!isNeighbor) {
                    // Calculate velocity differences
                    CUDA_REAL dvx = sh_vx[tid] - pi_vx;
                    CUDA_REAL dvy = sh_vy[tid] - pi_vy;
                    CUDA_REAL dvz = sh_vz[tid] - pi_vz;

                    CUDA_REAL magnitude1 = dx * dvx + dy * dvy + dz * dvz;
                    CUDA_REAL inv_sqrt_m0 = rsqrt(magnitude0);
                    CUDA_REAL inv_m0 = inv_sqrt_m0 * inv_sqrt_m0; // or 1 / magnitude0

                    CUDA_REAL scale = sh_mass[tid] * inv_sqrt_m0 * inv_m0;
                    ax = scale * dx;
                    ay = scale * dy;
                    az = scale * dz;

                    CUDA_REAL common_factor = 3.0 * magnitude1 * inv_m0;
                    adotx = scale * (dvx - common_factor * dx);
                    adoty = scale * (dvy - common_factor * dy);
                    adotz = scale * (dvz - common_factor * dz);
                }

                // Store results in the diff array
                diff[idx] = ax;
                diff[idx + n * m] = ay;
                diff[idx + 2 * n * m] = az;
                diff[idx + 3 * n * m] = adotx;
                diff[idx + 4 * n * m] = adoty;
                diff[idx + 5 * n * m] = adotz;
            }
            __syncthreads(); // Ensure all threads are done before loading the next tile
        }
    }
}
#endif