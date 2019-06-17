#include<stdio.h>
#include<sys/time.h>
#include "DEFINE_TYPE.h"
// Texture references for CSR matrix 
texture<int,1> tex_colm;
texture<int2,1> tex_val;


// Scratchpad used by vector dot product for reduction

FLOAT_TYPE* scratchpad;    





// Kernel block and grid parameters - threads in a block and blocks in a grid
#define NUM_THREADS 128
#define NUM_BLOCKS 16

// Macros to simplify kernels 
#define THREAD_ID threadIdx.x+blockIdx.x*blockDim.x
#define THREAD_COUNT gridDim.x*blockDim.x

// Solver parameters - relative tolerance and maximum iterations
#define epsilon 1e-7
#define IMAX 40000

// For timing solver

FLOAT_TYPE utime () {
  struct timeval tv;

  gettimeofday (&tv, NULL);

  return (tv.tv_sec + FLOAT_TYPE (tv.tv_usec) * 1e-6);
}

// Creates a diagonal matrix stored in a vector pcmat, from the CSR matrix findrm, colm, val.
// n is the matrix size.

__global__ void create_jac(int n, int* findrm, int* colm, FLOAT_TYPE* val, FLOAT_TYPE* pcmat)
{
  for(int i=THREAD_ID; i<n; i+=THREAD_COUNT) 
    for(int k=findrm[i]-1; k<findrm[i+1]-1; k++) 
      if(colm[k]==i+1)
        pcmat[i] = 1.0/val[k];
}

// Multiplies diagonal matrix mat stored as a vector by the vector src, storing result in dest.
// n is the vector length.


__global__ void diag_spmv(int n, FLOAT_TYPE *mat, FLOAT_TYPE *src, FLOAT_TYPE *dest) 
{
  for (int i=THREAD_ID; i<n; i+=THREAD_COUNT)
    dest[i] = mat[i]*src[i];
}

// Sets the length-n vector vec to the zero vector.

__global__ void veczero(int n, FLOAT_TYPE* vec) 
{
  for(int i=THREAD_ID; i<n; i+=THREAD_COUNT)
    vec[i] = 0;
}

// Allows fetching FLOAT_TYPE values from texture memory, which only supports integers

static __device__ FLOAT_TYPE fetch_FLOAT_TYPE(texture<int2,1> val, int elem)
{
  int2 v = tex1Dfetch(val, elem);
  return __hiloint2double(v.y, v.x);
}

// Multiplies the CSR matrix in findrm, tex_colm, tex_val by src and stores the
// result in dest. n is the matrix size/vector length.

__global__ void csr_spmv(int n, FLOAT_TYPE* src, FLOAT_TYPE* dest, int *findrm)
{
  for (int row=THREAD_ID; row<n; row+=THREAD_COUNT) {
    dest[row] = 0;
    int a=findrm[row];
    int b=findrm[row+1];
    for (int k=a;k<b;k++)
      dest[row] += fetch_FLOAT_TYPE(tex_val,k-1)*src[tex1Dfetch(tex_colm,k-1)-1];
  }
}

// Computes the dot product of length-n vectors vec1 and vec2. This is reduced in tmp into a
// single value per thread block. The reduced value is stored in the array partial.

__global__ void vecdot_partial(int n, FLOAT_TYPE* vec1, FLOAT_TYPE* vec2, FLOAT_TYPE* partial)
{ 
  __shared__ FLOAT_TYPE tmp[NUM_THREADS];
  tmp[threadIdx.x] = 0;

  for (int i=THREAD_ID; i<n; i+=THREAD_COUNT)
    tmp[threadIdx.x] += vec1[i]*vec2[i];

  for (int i=blockDim.x/2;i>=1;i = i/2) {
    __syncthreads();
    if (threadIdx.x < i) 
      tmp[threadIdx.x] += tmp[i + threadIdx.x]; 
  }

  if (threadIdx.x == 0) 
    partial[blockIdx.x] = tmp[0];
}

// Reduces the output of the vecdot_partial kernel to a single value. The result is stored in result.

__global__ void vecdot_reduce(FLOAT_TYPE* partial, FLOAT_TYPE* result)
{
  __shared__ FLOAT_TYPE tmp[NUM_BLOCKS];

  if (threadIdx.x < NUM_BLOCKS) 
    tmp[threadIdx.x] = partial[threadIdx.x]; 
  else 
    tmp[threadIdx.x] = 0;

  for (int i=blockDim.x/2;i>=1;i = i/2) {
    __syncthreads();
    if (threadIdx.x < i) 
      tmp[threadIdx.x] += tmp[i + threadIdx.x]; 
  }

  if (threadIdx.x == 0) 
    *result = tmp[0];
}

// Divides num by den and stores the result in result. This is very wasteful of the GPU.

__global__ void scalardiv(FLOAT_TYPE* num, FLOAT_TYPE* den, FLOAT_TYPE* result) 
{
  if(threadIdx.x==0 && blockIdx.x==0)
    *result = (*num)/(*den);
}

// Computes r= a*x+y for n-length vectors x and y, and scalar a.

__global__ void axpy(int n, FLOAT_TYPE* a, FLOAT_TYPE* x, FLOAT_TYPE* y, FLOAT_TYPE* r) 
{
  for (int i=THREAD_ID; i<n; i+=THREAD_COUNT)
    r[i] = y[i] + (*a)*x[i];
}

// Computes y= y-a*x for n-length vectors x and y, and scalar a.

__global__ void ymax(int n, FLOAT_TYPE* a, FLOAT_TYPE* x, FLOAT_TYPE* y) 
{
  for (int i=THREAD_ID; i<n; i+=THREAD_COUNT)
    y[i] = y[i] - (*a)*x[i];
}

// Convenient function for performing a vector dot product and reduce all in one go.

void vecdot(int n, FLOAT_TYPE* vec1, FLOAT_TYPE* vec2, FLOAT_TYPE* result) 
{ 
  dim3 BlockDim(NUM_THREADS);
  dim3 GridDim(NUM_BLOCKS);

  vecdot_partial<<<GridDim,BlockDim>>>(n, vec1, vec2, scratchpad);
  vecdot_reduce<<<1,NUM_BLOCKS>>>(scratchpad, result);
}

// Sets dest=src for scalars on the GPU.

void scalarassign(FLOAT_TYPE* dest, FLOAT_TYPE* src)
{
  cudaMemcpy(dest, src, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToDevice); 
}

// Sets dest=src for n-length vectors on the GPU.

void vecassign(FLOAT_TYPE *dest, FLOAT_TYPE *src, int n) 
{
  cudaMemcpy(dest, src, sizeof(FLOAT_TYPE)*n, cudaMemcpyDeviceToDevice);
}


// Main solver entry point. Solves Ax = b for x.
// CSR matrix (A) stores the row and column pointer information like Fortran indices (i.e. begin at 1)
// findrm_p     : matrix row pointer
// colm_p       : matrix column pointer
// matrix_val_p : matrix values
// size*        : size of each vector
// b_p          : pointer to RHS vector
// x_p          : solutzion (x) is returned here
extern "C"
void gpucg_solve_(int* findrm_p, int *size_findrm, int* colm_p, int* size_colm, FLOAT_TYPE* matrix_val_p, int *matrix_val_size, 
                             FLOAT_TYPE* b_p, int* rhs_val_size, FLOAT_TYPE *x_p)
{
  // CSR Matrix on the GPU
  int *k_findrm, *k_colm;
  FLOAT_TYPE *k_val;
  // Vectors on the GPU
  FLOAT_TYPE *k_b, *k_x, *k_r, *k_d, *k_q, *k_s;
  // Diagonal matrix on the GPU (stored as a vector)
  FLOAT_TYPE* k_jac;
  // Scalars on the GPU
  FLOAT_TYPE  *k_alpha, *k_snew, *k_beta, *k_sold, *k_s0;

  // Scalars on the host
  FLOAT_TYPE t, s0, snew;
  int iterations = 0;

  // Begin timing
  t = -utime ();

  // Allocate space on the GPU for the CSR matrix and RHS vector, and copy from host to GPU
  cudaMalloc((void**)&k_findrm, sizeof(int)*(*size_findrm));
  cudaMemcpy(k_findrm, findrm_p, sizeof(int)*(*size_findrm), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&k_colm, sizeof(int)*(*size_colm));
  cudaMemcpy(k_colm, colm_p, sizeof(int)*(*size_colm), cudaMemcpyHostToDevice);
  cudaBindTexture(NULL, tex_colm, k_colm, sizeof(int)*(*size_colm));
  cudaMalloc((void**)&k_val, sizeof(FLOAT_TYPE)*(*matrix_val_size));
  cudaMemcpy(k_val, matrix_val_p, sizeof(FLOAT_TYPE)*(*matrix_val_size), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&k_b, sizeof(FLOAT_TYPE)*(*rhs_val_size));
  cudaMemcpy(k_b, b_p, sizeof(FLOAT_TYPE)*(*rhs_val_size), cudaMemcpyHostToDevice);

  // Allocate space for vectors on the GPU
  cudaMalloc((void**)&k_x, sizeof(FLOAT_TYPE)*(*rhs_val_size));
  cudaMalloc((void**)&k_r, sizeof(FLOAT_TYPE)*(*rhs_val_size));
  cudaMalloc((void**)&k_d, sizeof(FLOAT_TYPE)*(*rhs_val_size));
  cudaMalloc((void**)&k_q, sizeof(FLOAT_TYPE)*(*rhs_val_size));
  cudaMalloc((void**)&k_s, sizeof(FLOAT_TYPE)*(*rhs_val_size));
  cudaMalloc((void**)&k_jac, sizeof(FLOAT_TYPE)*(*rhs_val_size));
  cudaMalloc((void**)&k_alpha, sizeof(FLOAT_TYPE));
  cudaMalloc((void**)&scratchpad, sizeof(FLOAT_TYPE)*NUM_BLOCKS);
  cudaMalloc((void**)&k_snew, sizeof(FLOAT_TYPE)*NUM_BLOCKS);
  cudaMalloc((void**)&k_sold, sizeof(FLOAT_TYPE));
  cudaMalloc((void**)&k_beta, sizeof(FLOAT_TYPE));
  cudaMalloc((void**)&k_s0, sizeof(FLOAT_TYPE));

  // Dimensions of blocks and grid on the GPU
  dim3 BlockDim(NUM_THREADS);
  dim3 GridDim(NUM_BLOCKS);

  struct timeval t1, t2;    
  gettimeofday(&t1, 0);

  // Create diagonal preconditioning matrix (J = 1/diag(M)) 
  create_jac<<<1,BlockDim>>>(*rhs_val_size, k_findrm, k_colm, k_val, k_jac);

  // Bind the matrix to the texture cache - this was not done earlier as we modified the matrix
  cudaBindTexture(NULL, tex_val, k_val, sizeof(FLOAT_TYPE)*(*matrix_val_size)); 

  // Initialise result vector (x=0)
  veczero<<<1,BlockDim>>>(*rhs_val_size, k_x);

  // r=b-Ax (r=b since x=0), and d=M^(-1)r
  cudaMemcpy(k_r, k_b, sizeof(FLOAT_TYPE)*(*rhs_val_size), cudaMemcpyDeviceToDevice);
  diag_spmv<<<1,BlockDim>>>(*rhs_val_size, k_jac, k_r, k_d);

  // s0 = r.d
  vecdot(*rhs_val_size, k_r, k_d, k_s0);
  // snew = s0
  scalarassign(k_snew, k_s0);

  // Copy snew and s0 back to host so that host can evaluate stopping condition
  cudaMemcpy(&snew, k_snew, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
  cudaMemcpy(&s0, k_s0, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);

  // While i < imax and snew > epsilon^2*s0
  while (iterations < IMAX && snew > epsilon*epsilon*s0)
  {
    // q = Ad
    csr_spmv<<<GridDim,BlockDim>>>(*rhs_val_size, k_d, k_q, k_findrm);
    // alpha = snew/(d.q)
    vecdot(*rhs_val_size, k_d, k_q, k_alpha);
    scalardiv<<<1,1>>>(k_snew, k_alpha, k_alpha);
    // x = x + alpha*d
    axpy<<<GridDim,BlockDim>>>(*rhs_val_size, k_alpha, k_d, k_x, k_x);
    // r = r - alpha*q
    ymax<<<GridDim,BlockDim>>>(*rhs_val_size, k_alpha, k_q, k_r);
    // s = M^(-1)r
    diag_spmv<<<GridDim,BlockDim>>>(*rhs_val_size, k_jac, k_r, k_s);
    // sold = snew
    scalarassign(k_sold, k_snew);
    // snew = r.s
    vecdot(*rhs_val_size, k_r, k_s, k_snew);
    // beta = snew/sold
    scalardiv<<<1,1>>>(k_snew, k_sold, k_beta);
    // d = s + beta*d
    axpy<<<GridDim,BlockDim>>>(*rhs_val_size, k_beta, k_d, k_s, k_d);
    // Copy back snew so the host can evaluate the stopping condition
    cudaMemcpy(&snew, k_snew, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    // i = i+1
    iterations++;
  }

  gettimeofday(&t2, 0);
  FLOAT_TYPE time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;

  printf("Time of GPU:  %3.1f ms \n", time);

  // Copy result vector back from GPU
  cudaMemcpy(x_p, k_x, sizeof(FLOAT_TYPE)*(*rhs_val_size), cudaMemcpyDeviceToHost);

  // Clean up
  cudaUnbindTexture(tex_colm);
  cudaUnbindTexture(tex_val);
  cudaFree(k_findrm);
  cudaFree(k_colm);
  cudaFree(k_val);
  cudaFree(k_b);
  cudaFree(k_x);
  cudaFree(k_r);
  cudaFree(k_d);
  cudaFree(k_q);
  cudaFree(k_jac);
  cudaFree(k_alpha);
  cudaFree(k_snew);
  cudaFree(k_sold);
  cudaFree(k_beta);
  cudaFree(k_s0);
  cudaFree(scratchpad);

  // End timing - call cudaThreadSynchronize so we know all computation is finished before we stop the clock.
  cudaThreadSynchronize();
  t += utime ();

  // Interesting information
  printf("Iterations: %d \n", iterations);
  printf("CUDA error is: %s \n", cudaGetErrorString(cudaGetLastError()));
  printf("Solve time: %f seconds\n", t);
}

