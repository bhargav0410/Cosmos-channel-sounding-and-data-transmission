#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <thread>
#include <chrono>
#include <ctime>
#include <random>
#include "immintrin.h"
//#include "chan_sounder.h"

using namespace std::chrono;

#define PI 3.141592654

//Generates PN sequence of 1 and -1 using LFSR based on polynomial given
std::vector<int> pn_seq_gen(std::vector<int> polynomial, int length) {
    std::vector<int> start_state(polynomial.size()-1, 0), output(length);
    start_state[0] = 1;
    int temp_bit;
    std::cout << "Starting PN seq gen...\n";
    for (int i = 0; i < length; i++) {
        temp_bit = 0;
        for (int p = 0; p < polynomial.size() - 1; p++) {
            if (polynomial[p] > 0) {
                temp_bit += start_state[p];
                temp_bit = temp_bit % 2;
            } 
        }
        //std::cout << temp_bit << ",";
        for (int r = 0; r < start_state.size() - 1; r++) {
            start_state[r] = start_state[r+1];
        }
        start_state[start_state.size() - 1] = temp_bit;
        output[i] = 2*temp_bit - 1;
        std::cout << output[i] << ",";
    }
    std::cout << "\n";
    std::cout << "PN sequence generated...\n";
    return output;
}

//Performs linear cross-correlation with a PN sequence to find PN sequence peaks in the received sequence and returns the index at which PN sequence starts
int find_pn_seq(std::complex<float> *in_seq, int *pn_seq, int in_size, int pn_size, int thres) {
    if (in_size < pn_size) {
        std::cout << "Correlation not possible...\n";
    }
    //out_seq = (std::complex<float> *)malloc((size_t)((in_size - pn_size + 1) * sizeof(std::complex<float>)));
    float temp;
    for (int i = 0; i < in_size - pn_size + 1; i++) {
        temp = 0;
        for (int j = 0; j < pn_size; j++) {
            temp += std::abs(in_seq[i + j] * std::complex<float>((float)pn_seq[j], 0));
        }
        temp = temp/pn_size;
        if (temp > thres) {
            return i;
        }
    }
    return -1;
}

//FFT without any library
inline void win_fft(std::complex<float> *fft_in, std::complex<float> *fft_out, int fft_size, int flag) __attribute__((optimize("-O3")));
inline void win_fft(std::complex<float> *fft_in, std::complex<float> *fft_out, int fft_size, int flag) {
    if (fft_size == 1) {
        fft_out[0] = fft_in[0];
    } else {
        std::vector<std::complex<float>> in_even(fft_size/2), in_odd(fft_size/2), out_first_h(fft_size/2), out_second_h(fft_size/2);
        for (int i = 0; i < fft_size/2; i++) {
            in_even[i] = fft_in[i*2];
            in_odd[i] = fft_in[i*2 + 1];
        }
        win_fft(&in_even[0], &out_first_h[0], fft_size/2, flag);
        win_fft(&in_odd[0], &out_second_h[0], fft_size/2, flag);
        for (int i = 0; i < fft_size/2; i++) {
            //FT for current stage
            //std::complex<float> temp = out_first_h[i];
            fft_out[i] = out_first_h[i] + exp(std::complex<float>(0, (float)flag * 2.0 * PI * (float)i/(float)fft_size)) * out_second_h[i];
            fft_out[fft_size/2 + i] = out_first_h[i] - exp(std::complex<float>(0, (float)flag * 2.0 * PI * (float)i/(float)fft_size)) * out_second_h[i];
        }
    }
}

//Performs element by element multiplication of one complex vector and conjugate of another complex vector and stores answer in third vector
inline void mult_by_conj(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len) __attribute__((optimize("-O3")));
inline void mult_by_conj(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len) {
	for (int i = 0; i < len; i++) {
        out[i] = in_vec[i] * std::conj(conj_vec[i]);
    }
}

inline void mult_by_conj_avx(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len) __attribute__((optimize("-O3")));
inline void mult_by_conj_avx(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len) {
    float temp_real[8], temp_imag[8];
    for (int i = 0; i < len; i += 8) {
        if (i + 8 < len) {
            __m256 in_1_r = _mm256_set_ps(in_vec[i+7].real(), in_vec[i+6].real(), in_vec[i+5].real(), in_vec[i+4].real(), in_vec[i+3].real(), in_vec[i+2].real(), in_vec[i+1].real(), in_vec[i].real());
            __m256 in_1_i = _mm256_set_ps(in_vec[i+7].imag(), in_vec[i+6].imag(), in_vec[i+5].imag(), in_vec[i+4].imag(), in_vec[i+3].imag(), in_vec[i+2].imag(), in_vec[i+1].imag(), in_vec[i].imag());
            __m256 in_2_r = _mm256_set_ps(conj_vec[i+7].real(), conj_vec[i+6].real(), conj_vec[i+5].real(), conj_vec[i+4].real(), conj_vec[i+3].real(), conj_vec[i+2].real(), conj_vec[i+1].real(), conj_vec[i].real());
            __m256 in_2_i = _mm256_set_ps(conj_vec[i+7].imag(), conj_vec[i+6].imag(), conj_vec[i+5].imag(), conj_vec[i+4].imag(), conj_vec[i+3].imag(), conj_vec[i+2].imag(), conj_vec[i+1].imag(), conj_vec[i].imag());
            __m256 out_r = _mm256_add_ps(_mm256_mul_ps(in_1_r, in_2_r), _mm256_mul_ps(in_1_i, in_2_i));
            __m256 out_i = _mm256_sub_ps(_mm256_mul_ps(in_2_r, in_1_i), _mm256_mul_ps(in_2_i, in_1_r));
            _mm256_store_ps(&temp_real[0], out_r);
            _mm256_store_ps(&temp_imag[0], out_i);
            for (int j = 0; j < 8; j++) {
                out[i + j] = std::complex<float>(temp_real[j], temp_imag[j]);
            }
        } else {
            for (int j = i; j < len; j++) {
                out[j] = in_vec[j] * std::conj(conj_vec[j]);
            }
            
        }
    }
}

inline void lin_corr_fft(std::complex<float> *in_vec1, std::complex<float> *in_vec2, std::complex<float> *out, int len1, int len2, int threads) __attribute__((optimize("-O3")));
inline void lin_corr_fft(std::complex<float> *in_vec1, std::complex<float> *in_vec2, std::complex<float> *out, int len1, int len2, int threads) {
    int corr_len = len1 + len2 - 1;
    int new_len = corr_len;
    std::vector<std::complex<float>> temp_in1, temp_in2, temp_out;
    if ((float)ceil(log2((float)corr_len)) - (float)log2((float)corr_len) > 0.001) {
        new_len = (int)pow(2,(int)ceil(log2((float)corr_len)));
    }
    //std::cout << "New length: " << new_len << "\n";
    temp_in1.resize(new_len, 0);
    temp_in2.resize(new_len, 0);
    temp_out.resize(new_len, 0);
    memcpy((void *)&temp_in1[0], (const void *)&in_vec1[0], len1*sizeof(std::complex<float>));
    memcpy((void *)&temp_in2[0], (const void *)&in_vec2[0], len1*sizeof(std::complex<float>));
    //for (int i = 0; i < len1; i++) {
    //    temp_in1[i] = in_vec1[i];
   // }
   // for (int i = 0; i < len2; i++) {
   //     temp_in2[i] = in_vec2[i];
   // }
    win_fft(&temp_in1[0], &temp_in1[0], new_len, -1);
    win_fft(&temp_in2[0], &temp_in2[0], new_len, -1);
    //for (int i = 0; i < new_len; i++) {
    //    temp_in1[i] = temp_in1[i]/std::complex<float>(new_len,0);
    //    temp_in2[i] = temp_in2[i]/std::complex<float>(new_len,0);
    //}
    //std::cout << "FFT performed of inputs...\n";
    mult_by_conj(&temp_in1[0], &temp_in2[0], &temp_out[0], new_len);
    //std::cout << "Multiplication with conjugate done...\n";
    win_fft(&temp_out[0], &temp_out[0], new_len, 1);
    for (int i = 0; i < new_len; i++) {
        temp_out[i] = temp_out[i]/std::complex<float>(new_len,0);
    }
    for (int i = 0; i < new_len; i++) {
        temp_out[i] = temp_out[i]/std::complex<float>(len2,0);
    }
    //std::cout << "IFFT performed...\n";
    memcpy((void *)&out[0], (const void *)&temp_out[0], corr_len*sizeof(std::complex<float>));
    //for (int i = 0; i < corr_len; i++) {
    //    out[i] = temp_out[i];
    //}

}

void lin_corr(std::complex<float> *in_vec1, std::complex<float> *in_vec2, std::complex<float> *out, int len1, int len2, int threads) {
    int corr_len = len1 + len2 - 1;
    std::vector<std::complex<float>> temp_in1, temp_in2, temp_out;
    temp_in1.resize(len1 + 2*(len2-1));
    temp_in2.resize(len2);
    for (int i = 0; i < len1; i++) {
        temp_in1[i + len2 - 2] = in_vec1[i];
    }
    for (int i = 0; i < len2; i++) {
        temp_in2[i] = in_vec2[i];
    }
    for (int i = 0; i < len1; i++) {
        for (int j = 0; j < len2; j++) {
            out[i] += temp_in1[i + j]*temp_in2[j];
        }
        out[i] = out[i]/std::complex<float>(len2,0);  
    }
}

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

//Performs element by element division of complex vectors and stores answer in numerator
void divide_by_vec(std::complex<float> *numer, std::complex<float> *denom, std::complex<float> *out, int len) {
	for (int i = 0; i < len; i++) {
        out[i] = numer[i]/denom[i];
    }
}

//Performs element by element division of complex vectors and stores answer in numerator (uses AVX for SIMD)
void divide_by_vec_avx(std::complex<float> *numer, std::complex<float> *denom, std::complex<float> *out, int len) {
    float num_real[8], num_imag[8], den_real[8], den_imag[8];

	for (int i = 0; i < len; i += 8) {
        if (i + 8 < len) {
            /*for (int j = 0; j < 8; j++) {
                num_real[j] = numer[i+j].real();
                num_imag[j] = numer[i+j].imag();
                den_real[j] = denom[i+j].real();
                den_imag[j] = denom[i+j].imag();
            }*/
            printf("Loading input...\n");
            __m256 num_in_r = _mm256_set_ps(numer[i+7].real(), numer[i+6].real(), numer[i+5].real(), numer[i+4].real(), numer[i+3].real(), numer[i+2].real(), numer[i+1].real(), numer[i].real());
            __m256 num_in_i = _mm256_set_ps(numer[i+7].imag(), numer[i+6].imag(), numer[i+5].imag(), numer[i+4].imag(), numer[i+3].imag(), numer[i+2].imag(), numer[i+1].imag(), numer[i].imag());
            __m256 den_in_r1 = _mm256_set_ps(denom[i+7].real(), denom[i+6].real(), denom[i+5].real(), denom[i+4].real(), denom[i+3].real(), denom[i+2].real(), denom[i+1].real(), denom[i].real());
            __m256 den_in_i1 = _mm256_set_ps(denom[i+7].imag(), denom[i+6].imag(), denom[i+5].imag(), denom[i+4].imag(), denom[i+3].imag(), denom[i+2].imag(), denom[i+1].imag(), denom[i].imag());
            __m256 den_in_r2 = den_in_r1;
            __m256 den_in_i2 = den_in_i1;
            printf("Getting denom...\n");
            __m256 den_out = _mm256_add_ps(_mm256_mul_ps(den_in_r1, den_in_r2), _mm256_mul_ps(den_in_i1, den_in_i2));
            __m256 num_out_r = _mm256_div_ps(_mm256_add_ps(_mm256_mul_ps(num_in_r, den_in_r1), _mm256_mul_ps(num_in_i, den_in_i1)), den_out);
            __m256 num_out_i = _mm256_div_ps(_mm256_sub_ps(_mm256_mul_ps(num_in_i, den_in_r1), _mm256_mul_ps(num_in_r, den_in_i1)), den_out);
            _mm256_store_ps((float *)&num_real[0], num_out_r);
            _mm256_store_ps((float *)&num_imag[0], num_out_i);
            std::cout << "Division performed...\n";
            for (int j = 0; j < 8; j++) {
                out[i+j] = std::complex<float>(num_real[j], num_imag[j]);
            }
            std::cout << "Got output...\n";
        } else {
            for (int j = i; j < len; j++) {
                out[j] = numer[j]/denom[j];
            }
        }
    }
}


int main(int argc, char *argv[]) __attribute__((optimize("-O3")));
int main(int argc, char *argv[]) {

    duration<double> timediff;
    high_resolution_clock::time_point start, finish;
    double ccorr_time;
    int pn_len = 1023;
    std::vector<int> polynomial = {1,0,0,0,0,0,0,1,0,0,1};
    std::vector<int> output, output_2;
    std::vector<float> pn_test;

    wh_matrix wh_mat;
    wh_mat.create_walsh_mat(3);

    for (int i = 0; i < wh_mat.wh_mat.size(); i++) {
        printf("|| ");
        for (int j = 0; j < wh_mat.wh_mat[i].size(); j++) {
            printf("%d\t",wh_mat.wh_mat[i][j]);
        }
        printf(" ||\n");
    }

    // PN Sequence test
    output = pn_seq_gen(polynomial, pn_len);
    output_2.resize(output.size());
    /*
    for (int i = 0; i < output_2.size(); i++) {
        output_2[i] = output[(i + 10) % output.size()];
    }
    for (int i = 0; i < output.size(); i++) {
        output[i] = output[i] + output_2[i];
    }
    */

    pn_test.resize(pn_len);
    start = high_resolution_clock::now();
    for (int i = 0; i < pn_len; i++) {
        for (int j = 0; j < pn_len; j++) {
            //std::cout << (int)((j + i) % pn_len) << ",";
            pn_test[i] = pn_test[i] + (output[j] * output[(int)((j + i) % pn_len)]);
        }
        //std::cout << "\n";
        pn_test[i] = pn_test[i]/(float)pn_len;
    }
    finish = high_resolution_clock::now();

    std::cout << "Circular corr time: " << duration_cast<duration<double>>(finish - start).count() << "\n";

    std::vector<std::complex<float>> num_in(pn_len), den_in(pn_len), div_out1(pn_len), div_out2(pn_len);
    for (int i = 0; i < pn_len; i++) {
        num_in[i] = std::complex<float>(rand(), rand());
        den_in[i] = std::complex<float>(rand(), rand());
    }    
    start = high_resolution_clock::now();
    divide_by_vec(&num_in[0], &den_in[0], &div_out1[0], pn_len);
    finish = high_resolution_clock::now();


    std::cout << "Divide by vector time: " << duration_cast<duration<double>>(finish - start).count() << "\n";

    /*
    start = high_resolution_clock::now();
    divide_by_vec_avx(&num_in[0], &den_in[0], &div_out2[0], pn_len);
    finish = high_resolution_clock::now();

    std::cout << "Divide by vector AVX time: " << duration_cast<duration<double>>(finish - start).count() << "\n";
    std::complex<float> err;
    for (int i = 0; i < pn_len; i++) {
        err += (div_out2[i] - div_out1[i]);
    }
    err = err/std::complex<float>(pn_len,0);
    std::cout << "Error: " << err << "\n";

    for (int i = 0; i < pn_len; i++) {
        std::cout << " Div1: " << div_out1[i] << ",";
        std::cout << " Div2: " << div_out2[i] << ",";
    }
    std::cout << "\n";
    */

//    for (int i = 0; i < pn_len; i++) {
//        std::cout << pn_test[i] << ",";
//    }
//    std::cout.flush();
//    std::cout << std::endl << std::endl;

    //FFT test
    std::vector<std::complex<float>> fft_in(pn_len + 1), fft_out(pn_len + 1), ifft_out(pn_len + 1);
    for (int i = 0; i < pn_len; i++) {
        fft_in[i] = std::complex<float>(output[i], 0);
    }
    start = high_resolution_clock::now();
    win_fft(&fft_in[0], &fft_out[0], pn_len + 1, -1);
    win_fft(&fft_in[0], &fft_out[0], pn_len + 1, -1);
    for (int i = 0; i < pn_len + 1; i++) {
        fft_out[i] = fft_out[i]/std::complex<float>(pn_len + 1, 0);
        fft_out[i] = fft_out[i]*conj(fft_out[i]);
        //std::cout << fft_out[i] << ",";
    }
    win_fft(&fft_out[0], &ifft_out[0], pn_len + 1, 1);
    finish = high_resolution_clock::now();

    std::cout << "FFT corr time: " << duration_cast<duration<double>>(finish - start).count() << "\n";
    
    std::vector<std::complex<float>> corr_out(2*pn_len - 1);
    //Linear correlation test
    fft_in.resize(pn_len);
    start = high_resolution_clock::now();
    for (int i = 0; i < 100; i++) {
        lin_corr_fft(&fft_in[0], &fft_in[0], &corr_out[0], pn_len, pn_len, 1);
    }
    finish = high_resolution_clock::now();

    std::cout << "FFT Linear corr time: " << duration_cast<duration<double>>(finish - start).count()/(double)100 << "\n";
/*
    for (int i = 0; i < corr_out.size(); i++) {
        if ((float)std::abs(corr_out[i]) < 0.5) {
            std::cout << "0, ";
        } else {
            std::cout << corr_out[i] << ",";
        }
        
    }
    std::cout << "\n";
*/
    start = high_resolution_clock::now();
    lin_corr(&fft_in[0], &fft_in[0], &corr_out[0], pn_len, pn_len, 1);
    finish = high_resolution_clock::now();

    std::cout << "Linear corr time: " << duration_cast<duration<double>>(finish - start).count() << "\n";

/*
    for (int i = 0; i < corr_out.size(); i++) {
        if ((float)std::abs(corr_out[i]) < 0.5) {
            std::cout << "0, ";
        } else {
            std::cout << corr_out[i] << ",";
        }
        
    }
    std::cout << "\n";
    */

    //for (int i = 0; i < pn_len + 1; i++) {
    //    std::cout << ifft_out[i] << ",";
   // }
   // std::cout.flush();
   // std::cout << std::endl << std::endl;

    //win_fft(&fft_in[0], &fft_out[0], pn_len + 1, 1);

//    for (int i = 0; i < pn_len + 1; i++) {
//        fft_out[i] = fft_out[i]/std::complex<float>((pn_len + 1), 0);
//        std::cout << fft_out[i] << ",";
//    }

    return 0;
}