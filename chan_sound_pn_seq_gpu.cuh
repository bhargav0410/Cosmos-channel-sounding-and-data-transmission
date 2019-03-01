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

__global__ void correlate(float *tx_array, float *rx_array, float *out_array, int thres, int num_rx_ants, int pn_len, int num_pn_seqs) {

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
    void correlate_with_pn_seq(std::vector<std::vector<std::complex<float>>> &input_from_radio) {
        //Allocating thread blocks in GPU for processing
        if (threads_alloc == 0) {
            rx_length = input_from_radio[0].size();
            num_rx_ants = (int)input_from_radio.size();
            cudaDeviceProp devProp;
            cudaGetDeviceProperties(&devProp, 0);
            int threads_per_block  = devProp.maxThreadsPerBlock;
            //X dimension of threads
            thread_x = std::max(threads_per_block, pn_length);
            //Y dimension of threads is used if length of PN sequence is smaller than number of threads per block
            if ((pn_length + 1) < threads_per_block) {
                thread_y = threads_per_block/(pn_length + 1);
            } else {
                thread_y = 1;
            }

            //X dimension of blocks is used for correlation
            block_x = ceil((rx_length + pn_length - 1)/thread_y);

            //Y dimension is used if the PN length is longer than number of threads in single block
            block_y = num_rx_ants;

            //Z dimension of blocks is for number of rx antennas
            block_z = num_rx_ants;


        }
        //Allocating memory for receiving array in GPU
        if (rx_alloc == 0) {
            cudaMalloc((void **)rx_array, input_from_radio[0].size()*num_rx_ants*sizeof(*rx_array));
            rx_alloc = 1;
        }
        //Allocating memory to output array
        if (out_alloc == 0) {
            cudaMalloc((void **)out_array, (input_from_radio[0].size() + pn_length - 1)*num_rx_ants*num_tx_seqs*sizeof(*out_array));
        }
        //Copying received signal to GPU
        for (int i = 0; i < input_from_radio.size(); i++) {
            cudaMemcpy((void *)&rx_array[(int)input_from_radio[0].size()*i], (void *)&input_from_radio[i][0], input_from_radio[0].size()*sizeof(*tx_array), cudaMemcpyHostToDevice);
        }
    }

private:
    int pn_length, num_tx_seqs, num_rx_ants, thres, rx_length, rx_alloc = 0, tx_alloc = 0, out_alloc = 0, threads_alloc = 0;
    int thread_x, thread_y, block_x, block_y, block_z;
    cuFloatComplex *tx_array, *rx_array, *out_array;

};

#endif