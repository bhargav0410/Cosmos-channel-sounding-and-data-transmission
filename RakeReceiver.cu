#include "RakeReceiver.cuh"

RakeReceiver::RakeReceiver {}

RakeReceiver::~RakeReceiver {}

__global__ void findPeaks(cuFloatComplex* H, float thres1, int* numPeaks, cuFloatComplex* gains, int* delays, int numAnts, int length) {
	float thres = thres1;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	
	__shared__ int temp = 0;
	
	if (tid < length) {
		if (cuCabsf(H[tid]) < thres) {
			gains[temp] = H[tid];
			delays[temp] = threadIdx.x;
			atomicAdd(&temp,1);
			
		}
	}
	numPeaks[blockIdx.x] = temp;
}

__global__ void combineChanGains(float Hsqrd, int* numPeaks, cuFloatComplex* gains, int numAnts) {
	
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	extern __shared__ temp[];
	int tempID = threadIdx.x;
	
	for (int i = 1; i < numPeaks*numAnts; i = i*2) {
		if (threadIdx.x%(2*i) == 0) {
			temp[tempID] = cuCaddf(temp[tempID],temp[tempID+i]);
		}
	}
	
}

__global__ void multiplyWithChannelConj(cuFloatComplex* Yf, cuFloatComplex* Y, cuFloatComplex* H, int* numPeaks, int numAnts, int length) {
	
	int tid = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*blockDim.x*gridDim.x + blockIdx.z*gridDim.y*gridDim.x*blockDim.x; 
	
	if (blockIdx.x*blockDim.x <= (length - blockDim.x)) {
		if () {
			Yf[tid] = 
		}
	}
	
}