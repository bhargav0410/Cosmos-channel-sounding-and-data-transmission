#ifndef CHAN_SOUND_PN_SEQ_GPU_CUH
#define CHAN_SOUND_PN_SEQ_GPU_CUH

#include <string>
#include <cstdlib>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuComplex.h>
#include <fstream>
#include <complex>
#include <algorithm>
#include <cmath>

//GPU function to correlate with a PN sequence
__global__ void correlate(float *tx_array, float *rx_array, float *out_array, int thres, int num_rx_ants, int pn_len, int rx_arr_len_per_ant) {

    extern __shared__ cuFloatComplex temp[];

    //Copying received data into shared memory by multiplying with the input PN sequence
    if (threadIdx.x + blockIdx.x + blockDim.x*blockIdx.y <= rx_arr_len_per_ant) {
        temp[threadIdx.x + pn_len*threadIdx.y] = cuCmulf(tx_array[threadIdx.x + blockDim.x*blockIdx.y], rx_array[threadIdx.x + blockIdx.x + blockDim.x*blockIdx.y + rx_arr_len_per_ant*blockIdx.z]);
    }
    __syncthreads();

    //Array reduction of data present in shared memory
    for (int i = 1; i < blockDim.x; i = i*2) {
        if (threadIdx.x % (2*i) == 0 && (threadIdx.x+i < blockDim.x)) {
            temp[i] = cuCaddf(temp[threadIdx.x], temp[threadIdx.x + i]);
        }
        __syncthreads();
    }
    float *real_part = (float *)&out_array[blockIdx.x + rx_arr_len_per_ant*blockIdx.z];
    float *imag_part = real_part + 1;
    //Reducing over multiple blocks
    if (threadIdx.x == 0) {
        atomicAdd(real_part, cuCrealf(temp[threadIdx.x]));
        atomicAdd(imag_part, cuCimagf(temp[threadIdx.x]));
    }
}

class chan_sound_pn_seq_gpu {

public:
    chan_sound_pn_seq_gpu() {}

    chan_sound_pn_seq_gpu(int _pn_len, int _num_tx_seqs, int _num_rx_ants, int _thres) {
        pn_length = _pn_len;
        num_tx_seqs = _num_tx_seqs;
        num_rx_ants = _num_rx_ants;
        thres = _thres;
    }

    ~chan_sound_pn_seq_gpu() {
        if (tx_alloc > 0) {
            free(tx_array);
            tx_alloc = 0;
        }
        if (rx_alloc > 0) {
            free(rx_array);
            rx_alloc = 0;
        }
        if (out_alloc > 0) {
            free(out_array);
            out_alloc = 0;
        }
    }

    //Function to write input PN sequences into the GPU
    void get_pn_seq_from_file(std::string file, int num_seqs_in_file, int seq_len) {
        num_tx_seqs = num_seqs_in_file;
        pn_length = seq_len;
        //Allocating memory for Tx array consising of PN sequences in GPU
        cudaMalloc((void **)tx_array, pn_length*num_tx_seqs*sizeof(*tx_array));
        tx_alloc = 1;
        std::ifstream infile(file.c_str(), std::ifstream::binary);
        float *temp_arr;
        temp_arr = (float *)malloc(2*pn_length*num_tx_seqs*sizeof(*temp_arr));
        infile.read(temp_arr, 2*pn_length*num_tx_seqs*sizeof(*temp_arr));
        //Copying to GPU
        cudaMemcpy((void *)tx_array, (void *)temp_arr, pn_length*num_tx_seqs*sizeof(*tx_array), cudaMemcpyHostToDevice);
        //Freeing temporary array
        free(temp_arr);
    }

    //Function to correlate with input PN sequence matrix 
    void correlate_with_pn_seq(std::vector<std::vector<std::complex<float>>> &input_from_radio, std::vector<std::vector<std::vector<std::complex<float>>>> out_array_cpu) {
        
        //Allocating memory for receiving array in GPU
        if (rx_alloc == 0) {
            cudaMalloc((void **)rx_array, 2*(pn_length - 1) + input_from_radio[0].size()*num_rx_ants*sizeof(*rx_array));
            rx_length = 2*(pn_length - 1) + input_from_radio[0].size()*num_rx_ants*sizeof(*rx_array);
            cudaMemset((void *)rx_array, 0, (pn_length - 1)*sizeof(*rx_array));
            cudaMemset((void *)&rx_array[rx_length - pn_length + 1], 0, (pn_length - 1)*sizeof(*rx_array));
            rx_alloc = 1;
        }
        //Allocating memory to output array
        if (out_alloc == 0) {
            cudaMalloc((void **)out_array, (input_from_radio[0].size() + pn_length - 1)*num_rx_ants*num_tx_seqs*sizeof(*out_array));
        }
        
        //Allocating thread blocks in GPU for processing
        if (threads_alloc == 0) {
            //rx_length = input_from_radio[0].size();
            num_rx_ants = (int)input_from_radio.size();
            cudaDeviceProp devProp;
            cudaGetDeviceProperties(&devProp, 0);
            int threads_per_block  = devProp.maxThreadsPerBlock;
            //X dimension of threads
            thread_x = std::min(threads_per_block, pn_length);
            //Y dimension of threads is used if length of PN sequence is smaller than number of threads per block
            /*if ((pn_length + 1) < threads_per_block) {
                thread_y = threads_per_block/(pn_length + 1);
            } else {
                thread_y = 1;
            }*/

            //X dimension of blocks is used for correlation
            block_x = rx_length + pn_length - 1;

            //Y dimension is used if the PN length is longer than number of threads in single block
            if (pn_length > threads_per_block) {
                block_y = (pn_length + 1)/threads_per_block;
            } else {
                block_y = 1;
            }

            //Z dimension of blocks is for number of rx antennas
            block_z = num_rx_ants;
        }
        dim3 blockDims(thread_x, 1), gridDims(block_x, block_y, block_z);
        
        //Copying received signal to GPU
        for (int i = 0; i < input_from_radio.size(); i++) {
            cudaMemcpy((void *)&rx_array[(int)input_from_radio[0].size()*i + pn_length - 1], (void *)&input_from_radio[i][0], input_from_radio[0].size()*sizeof(*tx_array), cudaMemcpyHostToDevice);
        }

        //Making sure output array is initialized to zero
        cudaMemset((void *)out_array, 0, (input_from_radio[0].size() + pn_length - 1)*num_rx_ants*num_tx_seqs*sizeof(*out_array));

        cudaStream_t stream[1024];
        for (int tx_stream = 0; tx_stream < num_tx_seqs; tx_stream++) {
            cudaStreamCreate(&stream[tx_stream]);
            //Performing correlation
            correlate <<<gridDims, blockDims, std::min(threads_per_block, pn_length)*sizeof(*rx_array), stream[tx_stream]>>> (tx_array, rx_array, &out_array[(input_from_radio[0].size() + pn_length - 1)*num_rx_ants*tx_stream], thres, num_rx_ants, pn_length, rx_length);
        }


    }

    //Function to set parameters such as threshold, number of pn sequences, PN sequence length, etc.
    void set_params(int _pn_length = 0, int _num_tx_seqs = 0, int _num_rx_ants = 0, int _thres = 0, int _rx_length = 0) {
        if (_pn_length > 0) {
            pn_length = _pn_length;
        }
        if (_num_tx_seqs > 0) {
            num_tx_seqs = _num_tx_seqs;
        }
        if (_num_rx_ants > 0) {
            num_rx_ants = _num_rx_ants;
        }
        if (_thres > 0) {
            thres = _thres;
        }
        if (_rx_length > 0) {
            rx_length = _rx_length;
        }
    }

protected:
    int pn_length, num_tx_seqs, num_rx_ants, thres, rx_length, rx_alloc = 0, tx_alloc = 0, out_alloc = 0, threads_alloc = 0;
    int thread_x, thread_y, block_x, block_y, block_z;
    cuFloatComplex *tx_array, *rx_array, *out_array;

};

#endif