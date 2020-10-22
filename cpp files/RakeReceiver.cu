#include "RakeReceiver.cuh"

RakeReceiver::RakeReceiver {}

RakeReceiver::~RakeReceiver {}

__global__ void findPeaks(cuFloatComplex* H, float thres1, int* numPeaks, cuFloatComplex* gains, int* delays, int numAnts, int length) {
	float thres = thres1;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	
	__shared__ int temp = 0;
	
	if (tid < length) {
		if (cuCabsf(H[tid]) > thres) {
			gains[temp] = H[tid];
			delays[temp] = threadIdx.x;
			atomicAdd(&temp,1);
			
		}
	}
	numPeaks[blockIdx.x] = temp;
}

__global__ void findMaxNumPeaks(int* maxVal, int* numPeaks, int numAnts) {
	
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	int tempID = threadIdx.x;
	extern __shared__ int temp[];
	
	if (tid < numAnts) {
		temp[tempID] = numPeaks[tid];
	}
	
	for (int i = 1; i < blockDim.x; i = i*2) {
		if (threadIdx.x%(2*i) == 0 and (threadIdx.x + i) < blockDim.x) {
			if (temp[tempID+i] > temp[tempID]) {
				temp[tempID] = temp[tempID+i];
			}
		}
	}
	if (tempID == 0) {
		maxVal[blockIdx.x] = temp[tempID];
	}
}

__global__ void combineChannelGainsOfPeaks(float* Hsqrd, int* numVals1, cuFloatComplex* gains) {
	int numVals = *numVals1;
	
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	extern __shared__ float temp[];
	int tempID = threadIdx.x;
	
	if (tempID < numVals) {
		temp[tempID] = cuCabsf(gains[tid])*cuCabsf(gains[tid]);
	} else {
		temp[tempID] = 0;
	}
	
	for (int i = 1; i < blockDim.x; i = i*2) {
		if (threadIdx.x%(2*i) == 0 and (threadIdx.x + i) < blockDim.x) {
			temp[tempID] += temp[tempID+i];
		}
	}
	__syncthreads();
	if (tempID == 0) {
		Hsqrd[blockIdx.x] = temp[tempID];
	}
}

__global__ void combineChannelGainsOfAntennas(float* Hsqrd, int* numVals1) {
	int numVals = *numVals1;
	extern __shared__ float temp[];
	int tempID = threadIdx.x;
	
	if (tempID < numVals) {
		temp[tempID] = Hsqrd[tempID];
	} else {
		temp[tempID] = 0;
	}
	
	for (int i = 1; i < blockDim.x; i = i*2) {
		if (threadIdx.x%(2*i) == 0 and (threadIdx.x + i) < blockDim.x) {
			temp[tempID] += temp[tempID+i];
		}
	}
	__syncthreads();
	if (tempID == 0) {
		Hsqrd[tempID] = temp[tempID];
	}	
}

__global__ void multiplyWithChannelConj(cuFloatComplex* Yf, cuFloatComplex* Y, cuFloatComplex* gains, int* delays, int* numPeaks, int numAnts, int length) {
	
	int tid = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*blockDim.x*gridDim.x + blockIdx.z*gridDim.y*gridDim.x*blockDim.x;
	int chanID = blockIdx.y + blockIdx.z*gridDim.y;
	int recID = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.z*gridDim.x*blockDim.x;
	
	if (blockIdx.x*blockDim.x <= (length - blockDim.x)) {
		if (blockIdx.y <= numPeaks[blockIdx.z]) {
			Yf[tid] = cuCmulf(gains[chanID],Y[recID]);
		}
	}
	
}