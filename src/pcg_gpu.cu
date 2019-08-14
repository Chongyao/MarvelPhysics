#include "pcg_gpu.h"
#include <iostream>
#include "cuda_runtime.h"
// #include "device_functions.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <fstream>
#include <assert.h>
#include <cmath>
#include <vector>
#include "cublas_v2.h"

// #include <Eigen/Dense>

#define N_MAX 2000
#define NNZ_MAX 13000

__constant__ int IA[N_MAX+1];
__constant__ int JA[NNZ_MAX];

void readIAandJA(const int size_Matrix,const int size_nozeronumber,const int *IAtemp, const int *JAtemp)
{
	cudaMalloc((void**)&IA, sizeof(int)*(size_Matrix + 1));
	cudaMalloc((void**)&JA, sizeof(int)*size_nozeronumber);
	cudaMemcpyToSymbol(IA, IAtemp, sizeof(int)*(size_Matrix + 1));
	cudaMemcpyToSymbol(JA, JAtemp, sizeof(int)*size_nozeronumber);
}

template<typename DOUBLE>
__global__ void initialvalue(int N, DOUBLE *A, DOUBLE *B, DOUBLE *Minverse, DOUBLE *r, DOUBLE *z, DOUBLE *p)
{
	int blockId = blockIdx.y*gridDim.x + blockIdx.x;
	int tid = blockId * (blockDim.x *blockDim.y) + (threadIdx.y*blockDim.x) + threadIdx.x;
	while (tid < N)
	{
		int jtmp = IA[tid + 1] - IA[tid];
		for (int j = 0; j < jtmp; j++)
		{
			if (JA[j + IA[tid]] == tid)
			{
				Minverse[tid] = 1.0 / A[j + IA[tid]];
			}
		}
		r[tid] = B[tid];
		z[tid] = Minverse[tid] * r[tid];
		p[tid] = z[tid];
		tid += (gridDim.x*blockDim.x)*(gridDim.y*blockDim.y);
	}
}

template<typename DOUBLE>
__global__ void VectorAMUtiplyP(int N, DOUBLE *A, DOUBLE *p, DOUBLE *ap)
{
	int blockId = blockIdx.y*gridDim.x + blockIdx.x;
	int tid = blockId * (blockDim.x *blockDim.y) + (threadIdx.y*blockDim.x) + threadIdx.x;
	while (tid < N)
	{
		DOUBLE temp = 0;
		int jtemp;
		jtemp = IA[tid + 1] - IA[tid];
		for (int j = 0; j < jtemp; j++)
		{
			temp += A[j + IA[tid]] * p[JA[j + IA[tid]]];
		}
		ap[tid] = temp;
		tid += (gridDim.x*blockDim.x)*(gridDim.y*blockDim.y);
	}
}

template<typename DOUBLE>
__global__ void inerate_ak(DOUBLE *zr, DOUBLE *pap, DOUBLE *ak)
{
	if (threadIdx.x == 0 && blockIdx.x == 0)
	{
		*ak = (*zr) / (*pap);
	}
}

template<typename DOUBLE>
__global__ void inerate_x(int N, DOUBLE *p, DOUBLE *ak, DOUBLE *x)
{
	int blockId = blockIdx.y*gridDim.x + blockIdx.x;
	int tid = blockId * (blockDim.x *blockDim.y) + (threadIdx.y*blockDim.x) + threadIdx.x;
	while (tid < N)
	{
		x[tid] = x[tid] + (*ak) * p[tid];
		tid += (gridDim.x*blockDim.x)*(gridDim.y*blockDim.y);
	}
}

template<typename DOUBLE>
__global__ void inerate_r(int N, DOUBLE *ak, DOUBLE *ap, DOUBLE *r)
{
	int blockId = blockIdx.y*gridDim.x + blockIdx.x;
	int tid = blockId * (blockDim.x *blockDim.y) + (threadIdx.y*blockDim.x) + threadIdx.x;
	while (tid < N)
	{
		r[tid] = r[tid] - (*ak)*ap[tid];
		tid += (gridDim.x*blockDim.x)*(gridDim.y*blockDim.y);
	}
}

template<typename DOUBLE>
__global__ void inerate_z(int N, DOUBLE *Minverse, DOUBLE *r, DOUBLE *z)
{
	int blockId = blockIdx.y*gridDim.x + blockIdx.x;
	int tid = blockId * (blockDim.x *blockDim.y) + (threadIdx.y*blockDim.x) + threadIdx.x;
	while (tid < N)
	{
		z[tid] = Minverse[tid] * r[tid];
		tid += (gridDim.x*blockDim.x)*(gridDim.y*blockDim.y);
	}
}

template<typename DOUBLE>
__global__ void inerate_p(int N, DOUBLE *zrnew, DOUBLE *zr, DOUBLE *z, DOUBLE *p)
{
	int blockId = blockIdx.y*gridDim.x + blockIdx.x;
	int tid = blockId * (blockDim.x *blockDim.y) + (threadIdx.y*blockDim.x) + threadIdx.x;
	while (tid < N)
	{
		p[tid] = z[tid] + ((*zrnew) / (*zr))*p[tid];
		tid += (gridDim.x*blockDim.x)*(gridDim.y*blockDim.y);
	}
}

template<typename DOUBLE>
void PCG<DOUBLE>::function_pcg(const int Ntemp, const int NNZtemp, const DOUBLE *Atemp, const DOUBLE *Btemp, DOUBLE* x)
{

        DOUBLE *A;
        DOUBLE *B;
        int N;
	N = Ntemp;
	assert(N <= N_MAX);
	assert(NNZtemp <= NNZ_MAX);
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	dim3 block(32, 32);
	dim3 grid((N + block.x - 1) / block.x, (N + block.y - 1) / block.y);
	DOUBLE *dev_Minverse, *dev_r, *dev_z, *dev_p;
	DOUBLE *zr = new DOUBLE, *dev_zr;
	DOUBLE *dev_ap;
	DOUBLE *pap = new DOUBLE, *dev_pap;
	DOUBLE *ak = new DOUBLE, *dev_ak;
	DOUBLE *dev_x;
	DOUBLE *zrnew = new DOUBLE, *dev_zrnew;

	cudaMalloc((void**)&A, sizeof(DOUBLE)*NNZtemp);
	cudaMalloc((void**)&B, sizeof(DOUBLE)*N);
	cudaMalloc((void**)&dev_Minverse, sizeof(DOUBLE)*N);
	cudaMalloc((void**)&dev_r, sizeof(DOUBLE)*N);
	cudaMalloc((void**)&dev_z, sizeof(DOUBLE)*N);
	cudaMalloc((void**)&dev_p, sizeof(DOUBLE)*N);
	cudaMalloc((void**)&dev_zr, sizeof(DOUBLE));
	cudaMalloc((void**)&dev_ap, sizeof(DOUBLE)*N);
	cudaMalloc((void**)&dev_pap, sizeof(DOUBLE));
	cudaMalloc((void**)&dev_ak, sizeof(DOUBLE));
	cudaMalloc((void**)&dev_x, sizeof(DOUBLE)*N);
	cudaMalloc((void**)&dev_zrnew, sizeof(DOUBLE));
	
	cudaMemcpy(A, Atemp, sizeof(DOUBLE)*NNZtemp, cudaMemcpyHostToDevice);
	cudaMemcpy(B, Btemp, sizeof(DOUBLE)*N, cudaMemcpyHostToDevice);


	initialvalue << <grid, block >> > (N, A, B, dev_Minverse, dev_r, dev_z, dev_p);
	cublasStatus_t status;
	cublasHandle_t handle;
	status = cublasCreate(&handle);
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		std::cout << "CUBLAS对象实例化出错" << std::endl;
		getchar();
	}
	for (int i = 0; i < N; i++)
	{
		cublasDdot(handle, N, dev_z, 1, dev_r, 1, dev_zr);
		VectorAMUtiplyP << <grid, block >> > (N, A, dev_p, dev_ap);
		cublasDdot(handle, N, dev_ap, 1, dev_p, 1, dev_pap);
		inerate_ak << <grid, block >> > (dev_zr, dev_pap, dev_ak);
		inerate_x << <grid, block >> > (N, dev_p, dev_ak, dev_x);
		inerate_r << <grid, block >> > (N, dev_ak, dev_ap, dev_r);
		inerate_z << <grid, block >> > (N, dev_Minverse, dev_r, dev_z);
		cublasDdot(handle, N, dev_z, 1, dev_r, 1, dev_zrnew);
		cudaMemcpy(zrnew, dev_zrnew, sizeof(DOUBLE), cudaMemcpyDeviceToHost);
		if (sqrt(*zrnew) < 1.0e-8) break;
		inerate_p << <grid, block >> > (N, dev_zrnew, dev_zr, dev_z, dev_p);
	}
	cudaMemcpy(x, dev_x, sizeof(DOUBLE)*N, cudaMemcpyDeviceToHost);
        


	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float time;
	cudaEventElapsedTime(&time, start, stop);
	std::cout << time;
	cudaFree(A);
	cudaFree(B);
	cudaFree(dev_Minverse);
	cudaFree(dev_r);
	cudaFree(dev_z);
	cudaFree(dev_p);
	cudaFree(dev_zr);
	cudaFree(dev_ap);
	cudaFree(dev_pap);
	cudaFree(dev_ak);
	cudaFree(dev_x);
	cudaFree(dev_zrnew);
}

// template __global__ void initialvalue<double>;
// template __global__ void initialvalue<float>;
// template __global__ void VectorAMUtiplyP<double>;
// template __global__ void VectorAMUtiplyP<float>;
// template __global__ void inerate_ak<double>;
// template __global__ void inerate_ak<float>;
// template __global__ void inerate_x<double>;
// template __global__ void inerate_x<float>;
// template __global__ void inerate_r<double>;
// template __global__ void inerate_r<float>;
// template __global__ void inerate_z<double>;
// template __global__ void inerate_z<float>;
// template __global__ void inerate_p<double>;
// template __global__ void inerate_p<float>;

