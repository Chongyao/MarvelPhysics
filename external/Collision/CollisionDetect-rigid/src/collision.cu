//atomicAdd(XX, 1) -> atomicInc !!!!!!!
#define OUTPUT_TXT

// CUDA Runtime
#include <cuda_runtime.h>

#include <cuda_profiler_api.h>
#include <assert.h>

// Utilities and system includes
#include <helper_functions.h>  // helper for shared functions common to CUDA SDK samples
#include <helper_cuda.h>       // helper for CUDA error checking

#include "vec3.cuh"
#include "tools.cuh"
#include "box.cuh"
#include "tri3f.cuh"
#include "bvh.cuh"
#include "pair.cuh"
#include "tri-contact.cuh"

#include <math.h>
#include <stdarg.h>

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

#include <string>
#include<algorithm>

using namespace std;

typedef unsigned int uint;


typedef struct {
	uint  numFace, numVert;
	REAL3 *_dx, *_dx0;
	tri3f *_df;
	g_box *_dfBx;

	// init function
	void init()
	{
		numFace = 0;
		numVert = 0;
		_dx0 = _dx = NULL;
		_df = NULL;
		_dfBx = NULL;
	}

	void destroy()
	{
		if (_dx == NULL) return;

		checkCudaErrors(cudaFree(_dx));
		checkCudaErrors(cudaFree(_dx0));
		checkCudaErrors(cudaFree(_df));
		checkCudaErrors(cudaFree(_dfBx));
	}

	void computeWSdata(REAL thickness, bool ccd);
} g_mesh;

//=======================================================

cudaDeviceProp deviceProp;
extern void initPairsGPU();

void initGPU()
{
	int devID = 0;
	cudaGetDevice(&devID);
	checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

	initPairsGPU();
}

//=======================================================

g_mesh theCloth;
g_bvh* theBVH;
g_front theFront;
g_pair thePairs[2]; // potentially colliding pairs
g_pair retPairs; //results

//rky
int BVH_NUM;
//=======================================================
//rky
void init(int num){
	BVH_NUM = num;
	theBVH=new g_bvh[num];
	//cudaMalloc((void **)&theFront, num*sizeof(g_front));
	theFront.init();
}
//rky
void front_clear() {
	theFront.clear();
}

void initPairsGPU()
{
	//pairs[0].init(MAX_PAIR_NUM); // MAX_PAIR_NUM);
	thePairs[1].init(MAX_PAIR_NUM);
	retPairs.init(MAX_PAIR_NUM / 10);
}

void pushMesh2GPU(int  numFace, int numVert, void *faces, void *nodes)
{
	theCloth.init();

	theCloth.numFace = numFace;
	theCloth.numVert = numVert;

	cudaMalloc((void **)&theCloth._df, numFace*sizeof(tri3f));
	cudaMalloc((void **)&theCloth._dfBx, numFace*sizeof(g_box));
	cudaMalloc((void **)&theCloth._dx, numVert*sizeof(REAL3));
	cudaMalloc((void **)&theCloth._dx0, numVert*sizeof(REAL3));

	cudaMemcpy(theCloth._df, faces, sizeof(tri3f)*numFace, cudaMemcpyHostToDevice);
	cudaMemcpy(theCloth._dx, nodes, sizeof(REAL3)*numVert, cudaMemcpyHostToDevice);
	cudaMemcpy(theCloth._dx0, theCloth._dx, sizeof(REAL3)*numVert, cudaMemcpyDeviceToDevice);

	theCloth.computeWSdata(0, false);
}

void updateMesh2GPU(void *nodes,void *prenode)
{
	cudaMemcpy(theCloth._dx0, prenode, sizeof(REAL3)*theCloth.numVert, cudaMemcpyHostToDevice);
	cudaMemcpy(theCloth._dx, nodes, sizeof(REAL3)*theCloth.numVert, cudaMemcpyHostToDevice);
	
	theCloth.computeWSdata(0, true);
}

//=======================================================

void pushBVHIdx(int max_level, unsigned int *level_idx, int i)
{
	theBVH[i]._max_level = max_level;
	theBVH[i]._level_idx = new uint[max_level];
	memcpy(theBVH[i]._level_idx, level_idx, sizeof(uint)*max_level);
}

void pushBVH(unsigned int length, int *ids, int i)
{
	theBVH[i]._num = length;
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bvh, length*sizeof(int) * 2));
	checkCudaErrors(cudaMemcpy(theBVH[i]._bvh, ids, length*sizeof(int) * 2, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bxs, length*sizeof(g_box)));
	checkCudaErrors(cudaMemset(theBVH[i]._bxs, 0, length*sizeof(g_box)));
	theBVH[i].hBxs = NULL;

	//rky
	theBVH[i]._triBxs = theCloth._dfBx;
	//theBVH[i]._triBxs =NULL;
	theBVH[i]._triCones = NULL;
}

void pushBVHself(unsigned int length, int *ids, bool i)
{
	theBVH[i]._num = length;
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bvh, length * sizeof(int) * 2));
	checkCudaErrors(cudaMemcpy(theBVH[i]._bvh, ids, length * sizeof(int) * 2, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bxs, length * sizeof(g_box)));
	checkCudaErrors(cudaMemset(theBVH[i]._bxs, 0, length * sizeof(g_box)));
	theBVH[i].hBxs = NULL;

	//rky
	theBVH[i]._triBxs =i? theCloth._dfBx:NULL;
	//theBVH[i]._triBxs =NULL;
	theBVH[i]._triCones = NULL;
}

void pushBVHLeaf(unsigned int length, int *idf, int i)
{
	checkCudaErrors(cudaMalloc((void**)&theBVH[i]._bvh_leaf, length*sizeof(int)));
	checkCudaErrors(cudaMemcpy(theBVH[i]._bvh_leaf, idf, length*sizeof(int), cudaMemcpyHostToDevice));
}

//======================================================


void refitBVH_Serial(int bvh_id, int length)
{

	refit_serial_kernel << <1, 1, 0 >> >
		(theBVH[bvh_id]._bvh, theBVH[bvh_id]._bxs, theBVH[bvh_id]._triBxs,
		theBVH[bvh_id]._cones, theBVH[bvh_id]._triCones,
		length == 0 ? theBVH[bvh_id]._num : length);

	getLastCudaError("refit_serial_kernel");
	cudaThreadSynchronize();
}

void refitBVH_Parallel(int bvh_id, int st, int length)
{
	BLK_PAR(length);

	refit_kernel << < B, T >> >
		(theBVH[bvh_id]._bvh, theBVH[bvh_id]._bxs, theBVH[bvh_id]._triBxs,
		theBVH[bvh_id]._cones, theBVH[bvh_id]._triCones,
		st, length);

	getLastCudaError("refit_kernel");
	cudaThreadSynchronize();
}

//rky
void refitBVH(int bvh_id)
{
	// before refit, need to get _tri_boxes !!!!
	// copying !!!
	for (int i = theBVH[bvh_id]._max_level - 1; i >= 0; i--) {
		int st = theBVH[bvh_id]._level_idx[i];
		int ed = (i != theBVH[bvh_id]._max_level - 1) ?
			theBVH[bvh_id]._level_idx[i + 1] - 1 : theBVH[bvh_id]._num - 1;

		int length = ed - st + 1;
		if (i < 5) {
			refitBVH_Serial(bvh_id, length + st);
			break;
		}
		else
		{
			refitBVH_Parallel(bvh_id, st, length);
		}
	}
}


void refitBVH()
{
	// before refit, need to get _tri_boxes !!!!
	// copying !!!
	for (int k = 0; k < BVH_NUM; k++){
		int bvh_id = k;
		for (int i = theBVH[bvh_id]._max_level - 1; i >= 0; i--) {
			int st = theBVH[bvh_id]._level_idx[i];
			int ed = (i != theBVH[bvh_id]._max_level - 1) ?
				theBVH[bvh_id]._level_idx[i + 1] - 1 : theBVH[bvh_id]._num - 1;

			int length = ed - st + 1;
			if (i < 5) {
				refitBVH_Serial(bvh_id, length + st);
				break;
			}
			else
			{
				refitBVH_Parallel(bvh_id, st, length);
			}
		}
	}
}

//===============================================

void pushFront(int num, unsigned int *data,unsigned int *data_id)
{
	g_front *f = &theFront;

	//rky
	//f->init();
	f->push(num, (uint4 *)data,data_id);
}

/*
void pushFrontself(bool self,int num, unsigned int *data)
{
	g_front *f = &theFront[self];

	//rky
	f->init();
	f->push(num, (uint4 *)data);
}
*/
//===============================================
// show memory usage of GPU
void  reportMemory(char *tag)
{
	//return;

#ifdef OUTPUT_TXT
	size_t free_byte;
	size_t total_byte;
	cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

	if (cudaSuccess != cuda_status) {
		printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status));
		exit(1);
	}

	REAL free_db = (REAL)free_byte;
	REAL total_db = (REAL)total_byte;
	REAL used_db = total_db - free_db;
	printf("%s: GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
		tag, used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
#endif
}

//===============================================

#define STACK_SIZE 1000
#define EMPTY (nIdx == 0)

#define PUSH_PAIR(nd1, nd2)  {\
	nStack[nIdx].x = nd1;\
	nStack[nIdx].y = nd2;\
	nIdx++;\
}

#define POP_PAIR(nd1, nd2) {\
	nIdx--;\
	nd1 = nStack[nIdx].x;\
	nd2 = nStack[nIdx].y;\
}

#define NEXT(n1, n2) 	POP_PAIR(n1, n2)


inline __device__ void pushToFront(int a, int b, uint4 *front, uint *idx, uint ptr, uint BVHA, uint BVHB, uint *bvh_id)
{
	//	(*idx)++;
	if (*idx < MAX_FRONT_NUM)
	{
		uint offset = atomicAdd(idx, 1);
		front[offset] = make_uint4(a, b, 0, ptr);
		
		bvh_id[offset * 2] = BVHA;
		bvh_id[offset * 2 + 1] = BVHB;
	}
}

//rky
inline __device__ void pushToBVHID(uint BVHA, uint BVHB,uint *idx,uint *bvh_id)
{
	if (*idx < MAX_FRONT_NUM * 2)
	{
		uint offset = *idx;
		bvh_id[offset * 2] = BVHA;
		bvh_id[offset * 2 + 1] = BVHB;
	}
}

inline __device__ void deleteFront_BVHID(uint4* front, uint* bvh_id, int idx,uint* frontIdx)
{
	uint tail = *frontIdx;

	//uint tembuffer[MAX_FRONT_NUM];
	/*
	cudaMemcpy(tembuffer, front+idx+1, offset * sizeof(uint), cudaMemcpyDeviceToDevice);
	cudaMemcpy(front + idx, tembuffer, offset * sizeof(uint), cudaMemcpyDeviceToDevice);

	cudaMemcpy(tembuffer, bvh_id + idx*2 + 2, offset*2 * sizeof(uint), cudaMemcpyDeviceToDevice);
	cudaMemcpy(bvh_id+idx*2, tembuffer, offset*2 * sizeof(uint), cudaMemcpyDeviceToDevice);
	*/

	for (int st = idx; st < tail; st++)
	{
		front[st] = front[st + 1];
	}

	for (int st = idx * 2; st < tail * 2; st++)
	{
		bvh_id[st] = bvh_id[st + 2];
	}

	atomicAdd(frontIdx, -1);
}

inline __device__ void First_P_node(uint left, uint right, int*bvhA, int* bvhB, int A_nums, int B_nums, uint bvh_idA, uint bvh_idB,
	g_box* bxsA, g_box* bxsB, uint4* front, uint* bvh_id, uint* frontidx,int idx)
{
	int ptr = front[idx].w;
	if (overlaps(left, right, bxsA, bxsB))
	{
		//pushToFront(left, right, front, frontidx, ptr);
		//pushToBVHID(bvh_idA, bvh_idB, frontidx, bvh_id);
	}
	else
	{

	}
}

inline __device__ void Front_Shrink(uint left, uint right, int* bvhA, int*bvhB, int A_nums,int B_nums, uint bvh_idA, uint bvh_idB,
	g_box* bxsA, g_box* bxsB, uint4* front, uint* bvh_id, uint* frontidx,int idx,int ptr)
{
	uint length = *frontidx;

	uint p_left = left == 0 ? 0 : getParent(left, bvhA, A_nums);
	uint p_right = right == 0 ? 0 : getParent(right, bvhB, B_nums);
	uint tem_left = left;
	uint tem_right = right;
	/*
	while (1)
	{		
		if (overlaps(p_left, p_right, bxsA, bxsB))
		{
			if (tem_left == getLeftChild(p_left, bvhA) && tem_right == getLeftChild(p_right, bvhB))
			{
				pushToFront(p_left, p_right, front, frontidx, ptr, bvh_idA, bvh_idB, bvh_id);
				//pushToBVHID(bvh_idA, bvh_idB, frontidx, bvh_id);
			}
			break;
		}
		else
		{
			if (p_left == 0 && p_right == 0)
			{
				if (tem_left == getLeftChild(p_left, bvhA) && tem_right == getLeftChild(p_right, bvhB))
				{
					pushToFront(p_left, p_right, front, frontidx, ptr, bvh_idA, bvh_idB, bvh_id);
					//pushToBVHID(bvh_idA, bvh_idB, frontidx, bvh_id);
				}
				break;
			}

			tem_left = p_left;
			tem_right = p_right;

			p_left = p_left == 0 ? 0 : getParent(p_left, bvhA, A_nums);
			p_right = p_right == 0 ? 0 : getParent(p_right, bvhB, B_nums);
		}
	}
	*/

	if (!overlaps(p_left, p_right, bxsA, bxsB))
	{
		while (1)
		{
			//��������Ǹ��ڵ� �ټ���Ƿ�������
			if (tem_left != p_left&&tem_right != p_right)
			{
				if (tem_left != getLeftChild(p_left, bvhA) || tem_right != getLeftChild(p_right, bvhB))
				{
					break;
				}
			}

			//�����һ�����ڵ� ��һ���������� break
			if (tem_left == p_left&&tem_right != getLeftChild(p_right, bvhB))
			{
				break;
			}

			if (tem_left != getLeftChild(p_left, bvhA) &&tem_right ==p_right)
			{
				break;
			}


			if (overlaps(p_left, p_right, bxsA, bxsB))
			{
				pushToFront(p_left, p_right, front, frontidx, ptr, bvh_idA, bvh_idB, bvh_id);
				break;
			}

			if (p_left == 0 && p_right == 0)
			{
				pushToFront(p_left, p_right, front, frontidx, ptr, bvh_idA, bvh_idB, bvh_id);
				break;
			}

			tem_left = p_left;
			tem_right = p_right;

			p_left = p_left == 0 ? 0 : getParent(p_left, bvhA, A_nums);
			p_right = p_right == 0 ? 0 : getParent(p_right, bvhB, B_nums);
		}
		
		front[idx].z = 1;
	}

	
	//deleteFront_BVHID(front, bvh_id, idx, frontidx);
}

/*
inline __device__ void sproutingAdaptive(int left, int right,
	int *bvhA, g_box *bxsA, int *bvhB, g_box *bxsB,
	uint4 *front, uint *frontIdx,
	uint2 *pairs, uint *pairIdx, bool update, uint ptr)
{
	uint2 nStack[STACK_SIZE];
	uint nIdx = 0;

	for (int i = 0; i<4; i++)
	{
		if (isLeaf(left, bvhA) && isLeaf(right, bvhB)) {
			pushToFront(left, right, front, frontIdx, ptr);
		}
		else {
			if (!overlaps(left, right, bxsA, bxsB)) {
				pushToFront(left, right, front, frontIdx, ptr);
			}
			else {
				if (isLeaf(left, bvhA)) {
					PUSH_PAIR(left, getLeftChild(right, bvhB));
					PUSH_PAIR(left, getRightChild(right, bvhB));
				}
				else {
					PUSH_PAIR(getLeftChild(left, bvhA), right);
					PUSH_PAIR(getRightChild(left, bvhA), right);
				}
			}
		}

		if (EMPTY)
			return;

		NEXT(left, right);
	}

	while (!EMPTY) {
		NEXT(left, right);
		pushToFront(left, right, front, frontIdx, ptr);
	}
}

*/
inline __device__ void sprouting(int left, int right,
	int *bvhA, g_box *bxsA, int *bvhB, g_box *bxsB,
	uint4 *front, uint *frontIdx,
	int2 *pairs, uint *pairIdx, bool update, uint ptr,uint *bvh_id,uint bvh_idA,uint bvh_idB)
{
	uint2 nStack[STACK_SIZE];
	uint nIdx = 0;

	while (1)
	{
		if (isLeaf(left, bvhA) && isLeaf(right, bvhB)) {
			if (update)
			{
				pushToFront(left, right, front, frontIdx, ptr, bvh_idA, bvh_idB, bvh_id);
				//pushToBVHID(bvh_idA, bvh_idB, frontIdx, bvh_id);
			}
			if (overlaps(left, right, bxsA, bxsB))
				addPair(getTriID(left, bvhA), getTriID(right, bvhB), pairs, pairIdx);
		}
		else {
			if (!overlaps(left, right, bxsA, bxsB)) {
				if (update)
				{
					pushToFront(left, right, front, frontIdx, ptr, bvh_idA, bvh_idB, bvh_id);
					//pushToBVHID(bvh_idA, bvh_idB, frontIdx, bvh_id);
				}
			}
			else {
				if (isLeaf(left, bvhA)) {
					PUSH_PAIR(left, getLeftChild(right, bvhB));
					PUSH_PAIR(left, getRightChild(right, bvhB));
				}
				else {
					PUSH_PAIR(getLeftChild(left, bvhA), right);
					PUSH_PAIR(getRightChild(left, bvhA), right);
				}
			}
		}

		if (EMPTY)
		{
			//cutilSafeCall(cudaMemcpy(rkylevel, %nx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
			//*rkylevel = nx;
			return;
		}
		NEXT(left, right);
	}
}

__device__ void doPropogate(
	uint4 *front,uint *bvh_id,g_bvh *bvh, uint *frontIdx, int num,
	int2 *pairs, uint *pairIdx, bool update, tri3f *Atris, int idx, bool *flags)
{
	uint4 node = front[idx];

	uint _bvhid[2];
	_bvhid[0] = bvh_id[idx*2];
	_bvhid[1] = bvh_id[idx*2 + 1];

	if (_bvhid[0] == _bvhid[1])
	{
		//front[idx].z = 1;
		//return;
	}

	if (node.z != 0) {
#if defined(_DEBUG) || defined(OUTPUT_TXT)
		//atomicAdd(frontIdx + 1, 1);
#endif
		//return;
	}

#ifdef USE_NC
	if (flags != NULL && flags[node.w] == 0) {
#if defined(_DEBUG) || defined(OUTPUT_TXT)
		atomicAdd(frontIdx + 2, 1);
#endif
		return;
	}
#endif


	uint left = node.x;
	uint right = node.y;
	int *bvhA = bvh[_bvhid[0]]._bvh;
	g_box *bxsA = bvh[_bvhid[0]]._bxs;
	int numsA = bvh[_bvhid[0]]._num;

	int *bvhB = bvh[_bvhid[1]]._bvh;
	g_box *bxsB = bvh[_bvhid[1]]._bxs;	
	int numsB = bvh[_bvhid[1]]._num;

	if (isLeaf(left, bvhA) && isLeaf(right, bvhB)) {
		if (overlaps(left, right, bxsA, bxsB))
		{
			if (_bvhid[0] != _bvhid[1])
				addPair(getTriID(left, bvhA), getTriID(right, bvhB), pairs, pairIdx);
			else { // for self ccd, we need to remove adjacent triangles, they will be processed seperatedly with orphan set
				if (!covertex(getTriID(left, bvhA), getTriID(right, bvhB), Atris))
					addPair(getTriID(left, bvhA), getTriID(right, bvhB), pairs, pairIdx);
			}
			return;
		}
	}
	
	int ptr = node.w;

	if (!overlaps(left, right, bxsA, bxsB))
	{
		if (update) 
		{
			if (left == 0 && right == 0)
				return;

			//Front_Shrink(left, right, bvhA, bvhB, numsA, numsB, _bvhid[0], _bvhid[1], bxsA, bxsB, front, bvh_id, frontIdx, idx, ptr);
		}

		return;
	}

	//if (update)
		//front[idx].z = 1;
	
	if (isLeaf(left, bvhA)) {
		sprouting(left, getLeftChild(right, bvhB), bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr,bvh_id,_bvhid[0], _bvhid[1]);
		sprouting(left, getRightChild(right, bvhB), bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr, bvh_id, _bvhid[0], _bvhid[1]);
	}
	else {
		sprouting(getLeftChild(left, bvhA), right, bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr, bvh_id, _bvhid[0], _bvhid[1]);
		sprouting(getRightChild(left, bvhA), right, bvhA, bxsA, bvhB, bxsB, front, frontIdx, pairs, pairIdx, update, ptr, bvh_id, _bvhid[0], _bvhid[1]);
	}

	if (update)
	{
		//deleteFront_BVHID(front, bvh_id, idx, frontIdx);
		front[idx].z = 1;
	}
}


__global__ void kernelPropogate(uint4 *front, uint *bvh_id,g_bvh *bvh, uint *frontIdx, int num,
	int2 *pairs, uint *pairIdx, bool update, tri3f *Atris, int stride, bool *flags)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	for (int i = 0; i<stride; i++) {
		int j = idx*stride + i;
		if (j >= num)
			return;

		doPropogate(front,bvh_id, bvh,frontIdx, num,
		 pairs, pairIdx, update, Atris, j, flags);
		
	}
}


int g_front::propogate(bool &update, bool ccd)
{
	uint dummy[1];
	cutilSafeCall(cudaMemcpy(dummy, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
#ifdef OUTPUT_TXT
	printf("Before propogate, length = %d\n", dummy[0]);
#endif

#if defined(_DEBUG) || defined(OUTPUT_TXT)
	uint dummy2[5] = { 0, 0, 0, 0, 0 };
	cutilSafeCall(cudaMemcpy(_dIdx + 1, dummy2, 5 * sizeof(int), cudaMemcpyHostToDevice));
#endif

	if (dummy[0] != 0) {
		//rky
		//g_bvh *pb1 = &theBVH[id1];
		//g_bvh *pb2 = &theBVH[id2];
		tri3f *faces = theCloth._df;

		int stride = 4;
#ifdef FIX_BT_NUM
		BLK_PAR2(dummy[0], stride);
#else
		BLK_PAR3(dummy[0], stride, getBlkSize((void *)kernelPropogate));
#endif
		g_bvh* _thebvh;
		cutilSafeCall(cudaMalloc((void**)&_thebvh, BVH_NUM*sizeof(g_bvh)));
		cutilSafeCall(cudaMemcpy(_thebvh, theBVH, BVH_NUM*sizeof(g_bvh),cudaMemcpyHostToDevice));
					
		double timing_start = omp_get_wtime();
			
		
		//cutilSafeCall(cudaMemcpy(rkylevel, theBVH, BVH_NUM * sizeof(g_bvh), cudaMemcpyHostToDevice));
		//rky
		kernelPropogate << < B, T >> >
			(_dFront,bvh_id,_thebvh, _dIdx, dummy[0],
			thePairs[1]._dPairs, thePairs[1]._dIdx, update, faces, stride,  
				theBVH[0]._ctFlags
				//&flags
				);
		
		double timing_finish = omp_get_wtime();
		
		/*
		FILE* file = fopen("out.txt", "a");
		fprintf(file, "%f\n", timing_finish-timing_start);
		fclose(file);
		*/

		/*
		cutilSafeCall(cudaMemcpy(dummy, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
		uint4 bvhid[50000];
		cutilSafeCall(cudaMemcpy(bvhid, _dFront, dummy[0] * sizeof(uint4), cudaMemcpyDeviceToHost));
		uint bvh_ids[50000];
		cutilSafeCall(cudaMemcpy(bvh_ids, bvh_id, dummy[0] * 2 * sizeof(uint), cudaMemcpyDeviceToHost));

		FILE* file = fopen("out.txt", "a");

		for (int st = 0; st < dummy[0]; st++)
		{
			fprintf(file, "%d %d %d %d %d\n", bvhid[st].x, bvh_ids[st * 2], bvhid[st].y, bvh_ids[st * 2 + 1], bvhid[st].z);
		}

		fprintf(file, "before sum: %d\n", dummy[0]);
		fprintf(file, "\n");
		*/
		cutilSafeCall(cudaMemcpy(dummy, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
		
		uint4* tem_front = new uint4[dummy[0]];
		uint* tem_bvhid = new uint[dummy[0] * 2];
		cutilSafeCall(cudaMemcpy(tem_front, _dFront, dummy[0] * sizeof(uint4), cudaMemcpyDeviceToHost));
		cutilSafeCall(cudaMemcpy(tem_bvhid, bvh_id, dummy[0]*2 * sizeof(uint), cudaMemcpyDeviceToHost));

		uint4* delete_front = new uint4[dummy[0]];
		uint* delete_bvhid = new uint[dummy[0] * 2];

		/*
		for (int st = 0; st < dummy[0]; st++)
		{
			for (int sk = st+1; sk < dummy[0]; sk++)
			{
				if (tem_front[st].x == tem_front[sk].x&&tem_front[st].y == tem_front[sk].y)
				{
					if (tem_bvhid[st * 2] == tem_bvhid[sk * 2] && tem_bvhid[st * 2 + 1] == tem_bvhid[sk * 2 + 1])
					{
						tem_front[sk].z = 1;
					}
				}
			}
		}
		*/
		double tem_timing_start = omp_get_wtime();
	
		int delete_nums = 0;
		for (int st = 0; st <dummy[0]; st++)
		{
			if (tem_front[st].z == 0)
			{				
				delete_front[delete_nums] = tem_front[st];
				delete_bvhid[delete_nums * 2] = tem_bvhid[st * 2];
				delete_bvhid[delete_nums * 2+1] = tem_bvhid[st * 2+1];

				delete_nums++;
			}
		}
		double tem_timing_finish = omp_get_wtime();
		
		/*
		FILE* file2 = fopen("out_delete.txt", "a");
		fprintf(file2, "%f\n", tem_timing_finish - tem_timing_start);
		fclose(file2);
		*/

		dummy[0] = delete_nums;
		cutilSafeCall(cudaMemcpy(_dIdx, dummy, 1 * sizeof(uint), cudaMemcpyHostToDevice));
		cutilSafeCall(cudaMemcpy(_dFront,delete_front, dummy[0] * sizeof(uint4), cudaMemcpyHostToDevice));
		cutilSafeCall(cudaMemcpy(bvh_id, delete_bvhid, dummy[0] * 2 * sizeof(uint), cudaMemcpyHostToDevice));

		delete[] tem_front;
		delete[] tem_bvhid;
		delete[] delete_front;
		delete[] delete_bvhid;

		
		//FILE* file = fopen("out.txt", "a");

		/*
		for (int st = 0; st < dummy[0]; st++)
		{
			fprintf(file, "%d %d %d %d %d\n", tem_front[st].x, tem_bvhid[st * 2], tem_front[st].y, tem_bvhid[st * 2 + 1], tem_front[st].z);
		}
		*/

		//fprintf(file, "%d\n",dummy[0]*sizeof(uint4));

		//fprintf(file, "after sum: %d\n", dummy[0]);
		//fprintf(file, "\n");

		//fclose(file);
		
		//thePairs[self]._dPairs, thePairs[self]._dIdx, update, faces, stride, (self && !ccd) ? theBVH[1]._ctFlags : NULL);
		//uint* rky = 0;
		//cutilSafeCall(cudaMemcpy(&rky, rkylevel, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
		//printf("rkylevel: %d\n", rky);
		cudaThreadSynchronize();
		getLastCudaError("kernelPropogate");
	}

	//cutilSafeCall(cudaMemcpy(dummy, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
#ifdef OUTPUT_TXT
	printf("After propogate, length = %d\n", dummy[0]);
#endif

#if defined(_DEBUG) || defined(OUTPUT_TXT)
	cutilSafeCall(cudaMemcpy(dummy2, _dIdx + 1, 5 * sizeof(int), cudaMemcpyDeviceToHost));
	//printf("Invalid = %d, NC culled = %d\n", dummy2[0], dummy2[1]);
#endif

	if (update && dummy[0] > SAFE_FRONT_NUM) {
		//printf("Too long front, stop updating ...\n");
		//update = false;
	}

	if (dummy[0] > MAX_FRONT_NUM) {
		printf("Too long front, exiting ...\n");
		exit(0);
	}
	return dummy[0];
}

//===============================================

__global__ void
kernel_face_ws(tri3f *face, REAL3 *x, REAL3 *ox, g_box *bxs, bool ccd, REAL thickness, int num)
{
	LEN_CHK(num);

	int id0 = face[idx].id0();
	int id1 = face[idx].id1();
	int id2 = face[idx].id2();

	REAL3 ox0 = ox[id0];
	REAL3 ox1 = ox[id1];
	REAL3 ox2 = ox[id2];
	REAL3 x0 = x[id0];
	REAL3 x1 = x[id1];
	REAL3 x2 = x[id2];

	bxs[idx].set(ox0, ox1);
	bxs[idx].add(ox2);

	if (ccd) {
		bxs[idx].add(x0);
		bxs[idx].add(x1);
		bxs[idx].add(x2);
	}
	//else
	bxs[idx].enlarge(thickness);
}

void g_mesh::computeWSdata(REAL thickness, bool ccd)
{
	if (numFace == 0)
		return;

	{
		int num = numFace;
		BLK_PAR(num);
		kernel_face_ws << <B, T >> > (
			_df, _dx, _dx0, _dfBx, ccd, thickness, num);
		getLastCudaError("kernel_face_ws");
	}
}

//===============================================

__device__ REAL triProduct(REAL3 &a, REAL3 &b, REAL3 &c)
{
	return dot(cross(a, b), c);
}

__device__ REAL3 xvpos(REAL3 x, REAL3 v, REAL t)
{
	return x + v*t;
}

__device__ int sgn(REAL x) { return x<0 ? -1 : 1; }

__device__ int solve_quadratic(REAL a, REAL b, REAL c, REAL x[2]) {
	// http://en.wikipedia.org/wiki/Quadratic_formula#Floating_point_implementation
	REAL d = b*b - 4 * a*c;
	if (d < 0) {
		x[0] = -b / (2 * a);
		return 0;
	}
	REAL q = -(b + sgn(b)*sqrt(d)) / 2;
	int i = 0;
	if (abs(a) > 1e-12*abs(q))
		x[i++] = q / a;
	if (abs(q) > 1e-12*abs(c))
		x[i++] = c / q;
	if (i == 2 && x[0] > x[1])
		fswap(x[0], x[1]);
	return i;
}

__device__ REAL newtons_method(REAL a, REAL b, REAL c, REAL d, REAL x0,
	int init_dir) {
	if (init_dir != 0) {
		// quadratic approximation around x0, assuming y' = 0
		REAL y0 = d + x0*(c + x0*(b + x0*a)),
			ddy0 = 2 * b + x0*(6 * a);
		x0 += init_dir*sqrt(abs(2 * y0 / ddy0));
	}
	for (int iter = 0; iter < 100; iter++) {
		REAL y = d + x0*(c + x0*(b + x0*a));
		REAL dy = c + x0*(2 * b + x0 * 3 * a);
		if (dy == 0)
			return x0;
		REAL x1 = x0 - y / dy;
		if (abs(x0 - x1) < 1e-6)
			return x0;
		x0 = x1;
	}
	return x0;
}


// solves a x^3 + b x^2 + c x + d == 0
__device__ int solve_cubic(REAL a, REAL b, REAL c, REAL d, REAL x[3]) {
	REAL xc[2];
	int ncrit = solve_quadratic(3 * a, 2 * b, c, xc);
	if (ncrit == 0) {
		x[0] = newtons_method(a, b, c, d, xc[0], 0);
		return 1;
	}
	else if (ncrit == 1) {// cubic is actually quadratic
		return solve_quadratic(b, c, d, x);
	}
	else {
		REAL yc[2] = { d + xc[0] * (c + xc[0] * (b + xc[0] * a)),
			d + xc[1] * (c + xc[1] * (b + xc[1] * a)) };
		int i = 0;
		if (yc[0] * a >= 0)
			x[i++] = newtons_method(a, b, c, d, xc[0], -1);
		if (yc[0] * yc[1] <= 0) {
			int closer = abs(yc[0])<abs(yc[1]) ? 0 : 1;
			x[i++] = newtons_method(a, b, c, d, xc[closer], closer == 0 ? 1 : -1);
		}
		if (yc[1] * a <= 0)
			x[i++] = newtons_method(a, b, c, d, xc[1], 1);
		return i;
	}
}

__device__ bool collision_test(
	const REAL3 &x0, const REAL3 &x1, const REAL3 &x2, const REAL3 &x3,
	const REAL3 &v0, const REAL3 &v1, const REAL3 &v2, const REAL3 &v3, REAL &time)
{
	REAL a0 = stp(x1, x2, x3),
		a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3),
		a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3),
		a3 = stp(v1, v2, v3);

	if (a1 == 0 && a2 == 0 && a3 == 0)
		return false;


	REAL t[4];
	int nsol = solve_cubic(a3, a2, a1, a0, t);
	//t[nsol] = 1; // also check at end of timestep

	bool t_flag = false;
	for (int i = 0; i < nsol; i++) {
		if (t[i] < 0 || t[i] > 1)
			continue;

		t_flag = true;
		if (t[i] < time)
			time = t[i];
	}
	return t_flag;
}

__device__ void doImpactVF(
	REAL3 x0, REAL3 x1, REAL3 x2, REAL3 x3,
	REAL3 x00, REAL3 x10, REAL3 x20, REAL3 x30, REAL &time
)
{
	REAL3 p0 = x00;
	REAL3 p1 = x10 - x00;
	REAL3 p2 = x20 - x00;
	REAL3 p3 = x30 - x00;
	REAL3 v0 = x0 - x00;
	REAL3 v1 = x1 - x10 - v0;
	REAL3 v2 = x2 - x20 - v0;
	REAL3 v3 = x3 - x30 - v0;

	/*	if (iii == 69 &&
	vid == 162 && t.id0() == 4 &&
	t.id1() == 44 && t.id2() == 32)
	vid = 162;
	*/
#ifdef  USE_DNF_FILTER
	//bool ret1 = dnf_filter(x00, x10, x20, x30, x0, x1, x2, x3);
#endif
	bool ret = collision_test(p0, p1, p2, p3, v0, v1, v2, v3, time);
	if (ret) {

	}
}

__device__ void doImpactEE(
	REAL3 x0, REAL3 x1, REAL3 x2, REAL3 x3,
	REAL3 x00, REAL3 x10, REAL3 x20, REAL3 x30, REAL &time)
{

	REAL3 p0 = x00;
	REAL3 p1 = x10 - x00;
	REAL3 p2 = x20 - x00;
	REAL3 p3 = x30 - x00;
	REAL3 v0 = x0 - x00;
	REAL3 v1 = x1 - x10 - v0;
	REAL3 v2 = x2 - x20 - v0;
	REAL3 v3 = x3 - x30 - v0;

#ifdef  USE_DNF_FILTER
	//bool ret1 = dnf_filter(x10, x20, x30, x40, x1, x2, x3, x4);
#endif
	/*	if (e0.x == 41 && e0.y == 624 &&
	e1.x == 599 && e1.y == 383)
	e0.x = 41;
	*/
	bool ret = collision_test(p0, p1, p2, p3, v0, v1, v2, v3, time);
	if (ret) {
	}
}


__global__ void kernelGetCollisions(
	int2 *pairs, int num, 
	REAL3 *cx, REAL3 *cx0, tri3f *ctris, int2 *pairRets, uint *pairIdx,double *t,
	int stride)
{
	int idxx = blockDim.x * blockIdx.x + threadIdx.x;

	for (int i = 0; i<stride; i++) {

		int j = idxx*stride + i;
		if (j >= num)
			return;

		int idx = j;

		int2 pair = pairs[idx];
		int fid1 = pair.x;
		int fid2 = pair.y;

		tri3f t1 = ctris[fid1];
		tri3f t2 = ctris[fid2];

		uint tt1[3];
		uint tt2[3];
		tt1[0] = t1.id0();
		tt1[1] = t1.id1();
		tt1[2] = t1.id2();

		tt2[0] = t2.id0();
		tt2[1] = t2.id1();
		tt2[2] = t2.id2();

		bool iscovetex = false;
		for (int st = 0; st < 3; st++)
		{
			for (int ss = 0; ss < 3; ss++)
			{
				if (tt1[st] == tt2[ss])
				{
					iscovetex = true;
				}
	}
}

		if (iscovetex)
			continue;

#ifdef FOR_DEBUG
		bool find = false;
		if (fid1 == 369 && fid2 == 3564)
			find = true;
		if (fid2 == 369 && fid1 == 3564)
			find = true;
#endif

		REAL3 p[3];
		p[0] = cx[t1.id0()];
		p[1] = cx[t1.id1()];
		p[2] = cx[t1.id2()];

		REAL3 pp[3];
		pp[0] = cx0[t1.id0()];
		pp[1] = cx0[t1.id1()];
		pp[2] = cx0[t1.id2()];

		REAL3 q[3];
		q[0] = cx[t2.id0()];
		q[1] = cx[t2.id1()];
		q[2] = cx[t2.id2()];

		REAL3 qq[3];
		qq[0] = cx0[t2.id0()];
		qq[1] = cx0[t2.id1()];
		qq[2] = cx0[t2.id2()];

		REAL time = 2;

		///*
		//VF
		for (int st = 0; st < 3; st++)
		{
		doImpactVF(
		p[st], q[0], q[1], q[2],
		pp[st], qq[0], qq[1], qq[2],
		time
		);
		}

		//VF
		for (int st = 0; st < 3; st++)
		{
		doImpactVF(q[st], p[0], p[1], p[2],
		qq[st], pp[0], pp[1], pp[2],
		time
		);
		}

		//EE
		int idx0[2];
		int idx1[2];
		for (int st = 0; st < 3; st++)
		{
		idx0[0] = st;
		idx0[1] = (st + 1) % 3;
		for (int ss = 0; ss < 3; ss++)
		{
		idx1[0] = ss;
		idx1[1] = (ss + 1) % 3;
		doImpactEE(
		cx[idx0[0]], cx[idx0[1]], cx[idx1[0]], cx[idx1[1]],
		cx0[idx0[0]], cx0[idx0[1]], cx0[idx1[0]], cx0[idx1[1]],
		time
		);
		}
		}

		if (time >= 0 && time <= 1)
		{
			int index=addPair(fid1, fid2,pairRets, pairIdx);
			if (index != -1)
			{
				t[index] = time;
			}
		}
		//*/


		/*
		REAL3 p0 = cx[t1.id0()];
		REAL3 p1 = cx[t1.id1()];
		REAL3 p2 = cx[t1.id2()];
		REAL3 q0 = cx[t2.id0()];
		REAL3 q1 = cx[t2.id1()];
		REAL3 q2 = cx[t2.id2()];

		if (tri_contact(p0, p1, p2, q0, q1, q2))
			addPair(fid1, fid2, pairRets, pairIdx);
		*/
	}
}

//===============================================

int g_pair::getCollisions(bool self, g_pair &ret,double *time)
{
	int num = length();
	printf("pair = %d\n", num);
#ifdef OUTPUT_TXT
	//if (self)
		//printf("self pair = %d\n", num);
	//else
		//printf("inter-obj pair = %d\n", num);
#endif

	if (num == 0)
		return 0;

	ret.clear();
	
	int stride = 4;
#ifdef FIX_BT_NUM
	BLK_PAR3(num, stride, 32);
#else
	BLK_PAR3(num, stride, getBlkSize((void *)kernelGetCollisions));
#endif
	
	REAL *tem_time;
	cutilSafeCall(cudaMalloc((void**)&tem_time, MAX_PAIR_NUM * sizeof(REAL)));

	kernelGetCollisions << < B, T >> >(_dPairs, num,
		theCloth._dx, theCloth._dx0, theCloth._df, ret._dPairs, ret._dIdx, tem_time,stride);

	getLastCudaError("kernelGetCollisions");

	int len = ret.length();

	if (len > 0)
	{
		cutilSafeCall(cudaMemcpy(time, tem_time, sizeof(REAL)* len, cudaMemcpyDeviceToHost));
	}

	cudaFree(tem_time);

#ifdef OUTPUT_TXT
	//printf("collision num = %d\n", len);
#endif

	return len;
}
//===============================================

int getCollisionsGPU(int *rets,double *time)
{
	bool update = true;
	int len = 0;
	
	double timing_start = omp_get_wtime();
	TIMING_BEGIN
	
	thePairs[1].clear();

	refitBVH();

	theFront.propogate(update,false);
	cudaThreadSynchronize();
	

	len = thePairs[1].getCollisions(true, retPairs,time);
	cudaThreadSynchronize();

	TIMING_END("$$$get_collisions_gpu")
	double timing_finish = omp_get_wtime();

	if (len > 0) {
		cutilSafeCall(cudaMemcpy(rets, retPairs._dPairs, sizeof(uint)*2*len, cudaMemcpyDeviceToHost));
	}

	//FILE* file = fopen("time.txt", "a");

	//fprintf(file, "%f\n", timing_finish - timing_start);

	//fclose(file);
	//fprintf(file, "\n");
	return len;
}
