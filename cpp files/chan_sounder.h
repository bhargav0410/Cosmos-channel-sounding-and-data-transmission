#ifndef CHAN_SOUNDER_H
#define CHAN_SOUNDER_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
//#include <fftw3.h>
#include <vector>
#include <complex>
#include <cmath>
#include <thread>
#include <omp.h>
#include "muFFT/fft.h"
#include "muFFT/fft_internal.h"
//#include "../muFFT/fft.h"
//#include "../muFFT/fft_internal.h"
#include <immintrin.h>


#define PI 3.141592654

//Walsh matrix struct
struct wh_matrix {
    std::vector<std::vector<int>> wh_mat;

    void create_walsh_mat(int power) {
        wh_mat.resize((int)pow(2,power));
        for (int i = 0; i < wh_mat.size(); i++) {
            wh_mat[i].resize((int)pow(2,power));
        }
        wh_mat[0][0] = 1;
        for (int i = 0; i < power; i++) {
            for (int ii = (int)pow(2,i); ii < (int)pow(2,i+1); ii++) {
                for (int jj = (int)pow(2,i); jj < (int)pow(2,i+1); jj++) {
                    wh_mat[ii - (int)pow(2,i)][jj] = wh_mat[ii - (int)pow(2,i)][jj - (int)pow(2,i)];
                    wh_mat[jj][ii - (int)pow(2,i)] = wh_mat[ii - (int)pow(2,i)][jj];
                    wh_mat[ii][jj] = -1*wh_mat[ii - (int)pow(2,i)][jj - (int)pow(2,i)];
                }
            }
        }
    }
};

//Distributed the tasks amongst workers with as much fairness as possible
void load_balancing_mpi(int *, int *, int, int) __attribute__((optimize("-O3")));

//Generates PN sequence of 1 and -1 using LFSR based on polynomial given
std::vector<int> pn_seq_gen(std::vector<int>, int) __attribute__((optimize("-O3")));

//Generates a zadoff-chu sequence of specific length
std::vector<std::complex<float>> zadoff_chu_gen(int) __attribute__((optimize("-O3")));

//Performs linear cross-correlation with a PN sequence to find PN sequence peaks in the received sequence and returns the index at which PN sequence starts
int find_pn_seq(std::complex<float> *, int *, int, int, float, int) __attribute__((optimize("-O3")));

int find_pn_seq_avx(std::complex<float> *, int *, int, int, float) __attribute__((optimize("-O3")));

//FFT without any library
void win_fft(std::complex<float> *, std::complex<float> *, int, int) __attribute__((optimize("-O3")));

//FFT of one row
void single_thread_fft(std::complex<float> *, std::complex<float> *, int) __attribute__((optimize("-O3")));

//IFFT of one row
void single_thread_ifft(std::complex<float> *, std::complex<float> *, int) __attribute__((optimize("-O3")));

/*Averages multiple vectors into one vector*/
void vector_averaging(std::complex<float> *, std::complex<float> *, int, int) __attribute__((optimize("-O3")));

//Swaps the [0:FFTsize/2-1] and [-FFTsize/2:FFTsize-1] halves of OFDM symbols and stores in same vector
void swap_halves(std::complex<float> *, int) __attribute__((optimize("-O3")));

//Performs element by element division of complex vectors and stores answer in numerator (uses AVX for SIMD)
void divide_by_vec_avx(std::complex<float> *, std::complex<float> *, std::complex<float> *, int) __attribute__((optimize("-O3")));

//Performs element by element multiplication of one complex vector and conjugate of another complex vector and stores answer in third vector
void mult_by_conj(std::complex<float> *, std::complex<float> *, std::complex<float> *, int) __attribute__((optimize("-O3")));

void mult_by_conj_avx(std::complex<float> *, std::complex<float> *, std::complex<float> *, int) __attribute__((optimize("-O3")));

void mult_by_conj_avx512(std::complex<float> *, std::complex<float> *, std::complex<float> *, int) __attribute__((optimize("-O3")));

//Finding maximum absolute value within vector
float find_max_val(std::complex<float> *, int, int) __attribute__((optimize("-O3")));

//Correlates one vector with cyclic shifts of another vector and gives output in a third vector.
//Total number of cyclic shifts are given by num_cyclic_shifts value. This value should include the zero cyclic shift which is basically correlation with a non-roatated vector.
void circ_correlate(std::complex<float> *, std::complex<float> *, std::complex<float> *, int, int) __attribute__((optimize("-O3")));

void circ_corr_fft(std::complex<float> *, std::complex<float> *, std::complex<float> *, int) __attribute__((optimize("-O3")));

void lin_corr_fft(std::complex<float> *, std::complex<float> *, std::complex<float> *, int, int) __attribute__((optimize("-O3")));

//Takes input vector of length fft_size - 1 and creates OFDM symbol of length fft_size + prefix_size
void create_ofdm_symbol(std::complex<float> *, std::complex<float> *, int, int) __attribute__((optimize("-O3")));

/*
Creates a channel sounding frame of OFDM symbols for multiple transmit antennas given by num_tx_ants.
This function assumes both inputs are pointers to C++ STL vectors. 
Also, there is an offset added of 'offset' number of samples before the OFDM symbols.
*/
void create_ofdm_sounding_frame(std::vector<std::complex<float>> &, std::vector<std::vector<std::complex<float>>> &, int, int, int, int, int, int, int) __attribute__((optimize("-O3")));

void demod_ofdm_symbol(std::complex<float> *, std::complex<float> *, int, int) __attribute__((optimize("-O3")));

void sound_frame(std::vector<std::complex<float>> &, std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &, int, int, int, int, int, int) __attribute__((optimize("-O3")));

void create_pn_seq_frame(std::vector<std::vector<int>>, std::vector<std::vector<std::complex<float>>> &, int, int, int, float, float, int) __attribute__((optimize("-O3")));

void sound_pn_frame(std::vector<std::vector<int>>, std::vector<std::vector<std::complex<float>>> &, std::vector<std::vector<std::complex<float>>> &, int, int, float, float, int) __attribute__((optimize("-O3")));

#endif