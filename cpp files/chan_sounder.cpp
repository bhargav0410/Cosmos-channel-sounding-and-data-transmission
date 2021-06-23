#include "chan_sounder.h"

//Distributed the tasks amongst workers with as much fairness as possible
void load_balancing_mpi(int *size_of_proc_data, int *displ, int num_procs, int len) {
    int displ_of_proc = 0;
    for (int i = 0; i < num_procs; i++) {
        displ[i] = displ_of_proc;
        size_of_proc_data[i] = (int)floor((float)len/(float)num_procs);
        if (i < (int)len % num_procs) {
            size_of_proc_data[i] += 1;
        }
        displ_of_proc += size_of_proc_data[i];
    }
}

//Generates PN sequence of 1 and -1 using LFSR based on polynomial given
std::vector<int> pn_seq_gen(std::vector<int> polynomial, int length) {
    std::vector<int> start_state(polynomial.size()-1, 0), output(length);
    start_state[0] = 1;
    int temp_bit;
    //std::cout << "Starting PN seq gen...\n";
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
        //std::cout << output[i] << ",";
    }
    //std::cout << "\n";
    //std::cout << "PN sequence generated...\n";
    return output;
}

//Generates a zadoff-chu sequence of specific length
std::vector<std::complex<float>> zadoff_chu_gen(int length) {
    std::vector<std::complex<float>> out(length);

    for (int i = 0; i < length; i++) {
        out[i] = std::exp(std::complex<float>(0,(float)(PI*i*(i + (length % 2)))/(float)length));
    }
    return out;
}

//Performs linear cross-correlation with a PN sequence to find PN sequence peaks in the received sequence and returns the index at which PN sequence starts
int find_pn_seq(std::complex<float> *in_seq, int *pn_seq, int in_size, int pn_size, float thres, int threads) {
    if (in_size < pn_size) {
        std::cout << "Correlation not possible...\n";
    }
    //out_seq = (std::complex<float> *)malloc((size_t)((in_size - pn_size + 1) * sizeof(std::complex<float>)));
    std::complex<float> temp;
    for (int i = 0; i < in_size - pn_size + 1; i++) {
        temp = 0;
        #pragma omp parallel num_threads(threads)
        {
            std::complex<float> local_temp = 0;
            #pragma omp for
            for (int j = 0; j < pn_size; j++) {
                local_temp += (in_seq[i + j] * std::complex<float>((float)pn_seq[j], 0));
            }
            //std::cout << "Local sum: " << local_temp << "\n";
            #pragma omp critical
            temp += local_temp;
        }
        //std::cout << temp << ",";
        if ((float)std::abs(temp)/(float)pn_size > thres) {
            return i;
        }
    }
    //std::cout << "\n";
    return -1;
}

int find_pn_seq_avx(std::complex<float> *in_seq, int *pn_seq, int in_size, int pn_size, float thres) {
    if (in_size < pn_size) {
        std::cout << "Correlation not possible...\n";
    }
    std::complex<float> temp;
    float in_real[8], in_imag[8], pn_real[8], temp_out[8];
    //in_real = (float *)malloc((size_t)(8 * sizeof(float)));
    //in_imag = (float *)malloc((size_t)(8 * sizeof(float)));
    //pn_real = (float *)malloc((size_t)(8 * sizeof(float)));
    //temp_out = (float *)malloc((size_t)(8 * sizeof(float)));
    int start = 0;
    for (int i = 0; i < in_size; i++) {
        if (std::abs(in_seq[i]) >= thres) {
            start = i;
            break;
        }
    }
    for (int i = start; i < in_size - pn_size + 1; i++) {
        temp = 0;
        for (int j = 0; j < pn_size; j += 8) {
            if (j + 8 < pn_size) {
                //std::cout << "Separating real and imag...\n";
                for (int jj = 0; jj < 8; jj++) {
                    in_real[jj] = in_seq[i + j + jj].real();
                    in_imag[jj] = in_seq[i + j + jj].imag();
                    pn_real[jj] = (float)pn_seq[j + jj];
                }
                //std::cout << "Loading in1...\n";
                __m256 in_r = _mm256_loadu_ps((const float *)&in_real[0]);
                //std::cout << "Loading in2...\n";
                __m256 in_i = _mm256_loadu_ps((const float *)&in_imag[0]);
                //std::cout << "Loading PN seq...\n";
                __m256 pn_r = _mm256_loadu_ps((const float *)&pn_real[0]);
                //std::cout << "Multiplying...\n";
                __m256 out_r = _mm256_mul_ps(in_r, pn_r);
                __m256 out_i = _mm256_mul_ps(in_i, pn_r);
                //std::cout << "Adding...\n";
                __m256 temp_add = _mm256_hadd_ps(out_r, out_i);
                temp_add = _mm256_hadd_ps(temp_add, temp_add);
                __m256 temp_perm = _mm256_permute2f128_ps(temp_add, temp_add, 0x01);
                temp_add = _mm256_add_ps(temp_add, temp_perm);
                _mm256_storeu_ps(&temp_out[0], temp_add);
                //std::cout << "Storing output...\n";
                temp += std::complex<float>(temp_out[0], temp_out[1]);
            } else {
                temp += (in_seq[i + j] * std::complex<float>((float)pn_seq[j], 0));
            }
        }
        //std::cout << temp << ",";
        if ((float)std::abs(temp)/(float)pn_size >= thres) {
            std::cout << "Corr val: " << (float)std::abs(temp)/(float)pn_size << "\n";
            return i;
        }
    }
    //std::cout << "\n";
    return -1;
}

//FFT without any library
void win_fft(std::complex<float> *fft_in, std::complex<float> *fft_out, int fft_size, int flag) {
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

//FFT of one row
void single_thread_fft(std::complex<float> *fft_in, std::complex<float> *fft_out, int fft_size) {
    //Setting up plan to execute
    //fftwf_plan plan;
    //plan = fftwf_plan_dft_1d(fft_size, (fftwf_complex *)fft_in, (fftwf_complex *)fft_out, FFTW_FORWARD, /*FFTW_MEASURE*/ FFTW_ESTIMATE);

    //Executing fft
    //fftwf_execute(plan);
    //Destroying plan
    //fftwf_destroy_plan(plan);

    mufft_plan_1d *muplan = mufft_create_plan_1d_c2c(fft_size, MUFFT_FORWARD, 1);
    mufft_execute_plan_1d(muplan, (cfloat *)fft_out, (cfloat *)fft_in);
    mufft_free_plan_1d(muplan);
}

//IFFT of one row
void single_thread_ifft(std::complex<float> *fft_in, std::complex<float> *fft_out, int fft_size) {
    //Setting up plan to execute
    //fftwf_plan plan;
    //plan = fftwf_plan_dft_1d(fft_size, (fftwf_complex *)fft_in, (fftwf_complex *)fft_out, FFTW_BACKWARD, /*FFTW_MEASURE*/ FFTW_ESTIMATE);

    //Executin ifft
    //fftwf_execute(plan);
    //Destroying plan
    //fftwf_destroy_plan(plan);

    mufft_plan_1d *muplan = mufft_create_plan_1d_c2c(fft_size, MUFFT_INVERSE, 1);
    mufft_execute_plan_1d(muplan, (cfloat *)fft_out, (cfloat *)fft_in);
    mufft_free_plan_1d(muplan);
}

/*Averages multiple vectors into one vector*/
void vector_averaging(std::complex<float> *input, std::complex<float> *output, int num_vectors, int vector_len) {

    //std::complex<float> *temp;
    //temp = (std::complex<float> *)malloc((size_t)(vector_len * sizeof(std::complex<float>)));
    //Output vector has first input vector
    for (int j = 0; j < vector_len; j++) {
        output[j] = input[j];
    }
    //Summing values of all vectors in one vector
    for (int i =  1; i < num_vectors; i++) {
        for (int j = 0; j < vector_len; j++) {
            output[j] = output[j] + input[j + i*vector_len];
        }
    }
    //Dividing by number of vectors
    for (int j = 0; j < vector_len; j++) {
        output[j] = output[j]/(float)num_vectors;
    }
}

//Swaps the [0:FFTsize/2-1] and [-FFTsize/2:FFTsize-1] halves of OFDM symbols and stores in same vector
void swap_halves(std::complex<float> *vec, int fft_size) {
    std::vector<std::complex<float>> temp((int)ceil((float)fft_size/(float)2));
    for (int i = 0; i < fft_size/2; i++) {
        temp[i] = vec[i];
        vec[i] = vec[i + fft_size/2];
        vec[i + fft_size/2] = temp[i];
    }
}

//Dividing 8-16 elements at same time
inline void divide_8_elems(std::complex<float> *in1, std::complex<float> *in2, std::complex<float> *out) {
    for (int i = 0; i < 8; i++) {
        out[i] = in1[i]/in2[i];
    }
}

inline void divide_16_elems(std::complex<float> *in1, std::complex<float> *in2, std::complex<float> *out) {
    for (int i = 0; i < 16; i++) {
        out[i] = in1[i]/in2[i];
    }
}

//Performs element by element division of complex vectors and stores answer in numerator
void divide_by_vec(std::complex<float> *numer, std::complex<float> *denom, std::complex<float> *out, int len) {

	for (int i = 0; i < len; i += 16) {
        if (i + 16 < len) {
            divide_16_elems(&numer[i], &denom[i], &out[i]);
        } else {
            for (int j = i, j < len; j++) {
                out[i] = numer[i]/denom[i];
            }
        }
    }
}

//Performs element by element division of complex vectors and stores answer in numerator (uses AVX for SIMD)
void divide_by_vec_avx(std::complex<float> *numer, std::complex<float> *denom, std::complex<float> *out, int len) {
    float num_real[8], num_imag[8], den_real[8], den_imag[8];

	for (int i = 0; i < len; i += 8) {
        if (i + 8 < len) {
            for (int j = 0; j < 8; j++) {
                num_real[j] = numer[i+j].real();
                num_imag[j] = numer[i+j].imag();
                den_real[j] = denom[i+j].real();
                den_imag[j] = denom[i+j].imag();
            }
            //printf("Loading input...\n");
            __m256 num_in_r = _mm256_load_ps((const float *)&num_real[0]);
            __m256 num_in_i = _mm256_load_ps((const float *)&num_imag[0]);
            __m256 den_in_r = _mm256_load_ps((const float *)&den_real[0]);
            __m256 den_in_i = _mm256_load_ps((const float *)&den_imag[0]);
            //printf("Getting denom...\n");
            __m256 den_out = _mm256_add_ps(_mm256_mul_ps(den_in_r, den_in_r), _mm256_mul_ps(den_in_i, den_in_i));
            __m256 num_out_r = _mm256_div_ps(_mm256_add_ps(_mm256_mul_ps(num_in_r, den_in_r), _mm256_mul_ps(num_in_i, den_in_i)), den_out);
            __m256 num_out_i = _mm256_div_ps(_mm256_sub_ps(_mm256_mul_ps(num_in_i, den_in_r), _mm256_mul_ps(num_in_r, den_in_i)), den_out);
            _mm256_store_ps((float *)&num_real[0], num_out_r);
            _mm256_store_ps((float *)&num_imag[0], num_out_i);
            for (int j = 0; j < 8; j++) {
                out[i+j] = std::complex<float>(num_real[j], num_imag[j]);
            }
        } else {
            for (int j = i; j < len; j++) {
                out[j] = numer[j]/denom[j];
            }
        }
    }
}

//Performs element by element multiplication of one complex vector and conjugate of another complex vector and stores answer in third vector
void mult_by_conj(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len) {
	for (int i = 0; i < len; i += 16) {
        if (i + 16 < len) {
            for (int j = 0; j < 16; j++) {
                out[i + j] = in_vec[i + j] * std::conj(conj_vec[i + j]);
            }
        } else {
            for (int j = i; j < len; j++) {
                out[j] = in_vec[j] * std::conj(conj_vec[j]);
            }
        }
    }
}

void mult_by_conj_avx(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len) {
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
            _mm512_store_ps(&temp_imag[0], out_i);
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

void mult_by_conj_avx512(std::complex<float> *in_vec, std::complex<float> *conj_vec, std::complex<float> *out, int len) {
    float temp_real[16], temp_imag[16];
    for (int i = 0; i < len; i += 16) {
        if (i + 16 < len) {
            __m512 in_1_r = _mm512_set_ps(in_vec[i+15].real(), in_vec[i+14].real(), in_vec[i+13].real(), in_vec[i+12].real(), in_vec[i+11].real(), in_vec[i+10].real(), in_vec[i+9].real(), in_vec[i+8].real(), in_vec[i+7].real(), in_vec[i+6].real(), in_vec[i+5].real(), in_vec[i+4].real(), in_vec[i+3].real(), in_vec[i+2].real(), in_vec[i+1].real(), in_vec[i].real());
            __m512 in_1_i = _mm512_set_ps(in_vec[i+15].imag(), in_vec[i+14].imag(), in_vec[i+13].imag(), in_vec[i+12].imag(), in_vec[i+11].imag(), in_vec[i+10].imag(), in_vec[i+9].imag(), in_vec[i+8].imag(), in_vec[i+7].imag(), in_vec[i+6].imag(), in_vec[i+5].imag(), in_vec[i+4].imag(), in_vec[i+3].imag(), in_vec[i+2].imag(), in_vec[i+1].imag(), in_vec[i].imag());
            __m512 in_2_r = _mm512_set_ps(conj_vec[i+15].real(), conj_vec[i+14].real(), conj_vec[i+13].real(), conj_vec[i+12].real(), conj_vec[i+11].real(), conj_vec[i+10].real(), conj_vec[i+9].real(), conj_vec[i+8].real(), conj_vec[i+7].real(), conj_vec[i+6].real(), conj_vec[i+5].real(), conj_vec[i+4].real(), conj_vec[i+3].real(), conj_vec[i+2].real(), conj_vec[i+1].real(), conj_vec[i].real());
            __m512 in_2_i = _mm512_set_ps(conj_vec[i+15].imag(), conj_vec[i+14].imag(), conj_vec[i+13].imag(), conj_vec[i+12].imag(), conj_vec[i+11].imag(), conj_vec[i+10].imag(), conj_vec[i+9].imag(), conj_vec[i+8].imag(), conj_vec[i+7].imag(), conj_vec[i+6].imag(), conj_vec[i+5].imag(), conj_vec[i+4].imag(), conj_vec[i+3].imag(), conj_vec[i+2].imag(), conj_vec[i+1].imag(), conj_vec[i].imag());
            __m512 out_r = _mm512_add_ps(_mm512_mul_ps(in_1_r, in_2_r), _mm512_mul_ps(in_1_i, in_2_i));
            __m512 out_i = _mm512_sub_ps(_mm512_mul_ps(in_2_r, in_1_i), _mm512_mul_ps(in_2_i, in_1_r));
            _mm512_store_ps(&temp_real[0], out_r);
            _mm512_store_ps(&temp_imag[0], out_i);
            for (int j = 0; j < 16; j++) {
                out[i + j] = std::complex<float>(temp_real[j], temp_imag[j]);
            }
        } else {
            for (int j = i; j < len; j++) {
                out[j] = in_vec[j] * std::conj(conj_vec[j]);
            }
            
        }
    }
}

//Finding maximum absolute value within vector
float find_max_val(std::complex<float> *in_vec, int len, int threads) {
    std::vector<float> abs_vec(len);
    float temp_ret;

    //Getting absolute value of complex number
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < len; i++) {
        abs_vec[i] = std::abs(in_vec[i]);
    }

    //Finding max value of absolute value vector
    for (int step = 1; step < len; step *= 2) {

        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < len; i += step*2) {
            if (i + step < len) {
                abs_vec[i] = std::max(abs_vec[i], abs_vec[i+step]);
            }
        }
    }
    return abs_vec[0];
}

//Correlates one vector with cyclic shifts of another vector and gives output in a third vector.
//Total number of cyclic shifts are given by num_cyclic_shifts value. This value should include the zero cyclic shift which is basically correlation with a non-roatated vector.
void circ_correlate(std::complex<float> *in_vec, std::complex<float> *cyclic_shift_vec, std::complex<float> *out, int num_cyclic_shifts, int len) {
    //Performing correlation
    for (int i = 0; i < num_cyclic_shifts; i++) {
        for (int j = 0; j < len; j++) {
            out[i] = out[i] + (in_vec[j] * cyclic_shift_vec[(int)((j + i) % num_cyclic_shifts)]);
        }
    }
}

void circ_corr_fft(std::complex<float> *in_vec, std::complex<float> *cyclic_shift_vec, std::complex<float> *out, int len) {
    int new_len = len;
    std::vector<std::complex<float>> temp_in1, temp_in2, temp_out;
    if ((float)ceil(log2((float)len)) - (float)log2((float)len) > 0.001) {
        new_len = (int)ceil(log2((float)len));
    }
    temp_in1.resize(new_len);
    temp_in2.resize(new_len);
    temp_out.resize(len);
    for (int i = 0; i < len; i++) {
        temp_in1[i] = in_vec[i];
        temp_in2[i] = cyclic_shift_vec[i];
    }
    single_thread_fft(&temp_in1[0], &temp_in1[0], new_len);
    single_thread_fft(&temp_in2[0], &temp_in2[0], new_len);
    mult_by_conj_avx(&temp_in1[0], &temp_in2[0], &temp_out[0], new_len);
    single_thread_ifft(&temp_out[0], &temp_out[0], new_len);
    for (int i = 0; i < len; i++) {
        out[i] = temp_out[i];
    }
}

void lin_corr_fft(std::complex<float> *in_vec1, std::complex<float> *in_vec2, std::complex<float> *out, int len1, int len2) {
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
    //std::cout << "Performing fft...\n";
    single_thread_fft(&temp_in1[0], &temp_in1[0], new_len);
    single_thread_fft(&temp_in2[0], &temp_in2[0], new_len);
    //std::cout << "Multiplying with conjugate...\n";
    mult_by_conj_avx(&temp_in1[0], &temp_in2[0], &temp_out[0], new_len);
    //std::cout << "Performing ifft...\n";
    single_thread_ifft(&temp_out[0], &temp_out[0], new_len);
    //std::cout << "Copying output...\n";
    memcpy((void *)&out[0], (const void *)&temp_out[0], corr_len*sizeof(std::complex<float>));

}

//Takes input vector of length fft_size - 1 and creates OFDM symbol of length fft_size + prefix_size
void create_ofdm_symbol(std::complex<float> *in_vec, std::complex<float> *out_vec, int fft_size, int prefix_size) {
    std::complex<float> *temp;
    temp = (std::complex<float> *)malloc((size_t)(fft_size * sizeof(std::complex<float>)));
    //Adding zero sub-carrier by copying input vector to temp vector
    memcpy((void *)&temp[0], (const void *)&in_vec[0], (int)floor((float)fft_size/(float)2)*sizeof(std::complex<float>));
    //temp[(int)floor((float)fft_size/(float)2)] = 0;
    memcpy((void *)&temp[(int)floor((float)fft_size/(float)2) + 1], (const void *)&in_vec[(int)floor((float)fft_size/(float)2)], ((int)floor((float)fft_size/(float)2)-1)*sizeof(std::complex<float>));
    //temp.insert(temp.begin() + fft_size/2, std::complex<float>(0,0));

    //Swapping halves
    swap_halves(&temp[0], fft_size);

    //Performing IFFT on input symbol
//    win_fft(&temp[0], &out_vec[prefix_size], fft_size, 1);
    single_thread_ifft(&temp[0], &out_vec[prefix_size], fft_size);

    //Adding cyclic prefix values
    for (int j = 0; j < prefix_size; j++) {
        out_vec[j] = out_vec[fft_size + j];
    }
    free(temp);
}

/*
Creates a channel sounding frame of OFDM symbols for multiple transmit antennas given by num_tx_ants.
This function assumes both inputs are pointers to C++ STL vectors. 
Also, there is an offset added of 'offset' number of samples before the OFDM symbols.
*/
void create_ofdm_sounding_frame(std::vector<std::complex<float>> &pilot_vec, std::vector<std::vector<std::complex<float>>> &out_vec, int fft_size, int prefix_size, int num_tx_ants, int pre_offset, int post_offset, int num_reps, int num_threads) {
    if (pilot_vec.size() != fft_size - 1) {
        std::cout << "Pilot vector size should be of the OFDM FFT size - 1";
        return;
    }
    if (out_vec.size() < num_tx_ants) {
        out_vec.resize(num_tx_ants);
    }
    for (int i = 0; i < num_tx_ants; i++) {
        out_vec[i].resize((int)(((fft_size + prefix_size)*num_tx_ants + pre_offset + post_offset)*num_reps));
    }
    std::vector<std::complex<float>> temp_out(fft_size + prefix_size);

    //std::cout << "Setting number of threads...\n";

    //omp_set_dynamic(0);
    //omp_set_num_threads(num_threads);

    //std::cout << "Using OpenMP for multi-threaded modulation...\n";
    //Initial modulation of pilot symbols
    create_ofdm_symbol(&pilot_vec[0], &temp_out[0], fft_size, prefix_size);
    //float max_val = find_max_val(&temp_out[0], fft_size + prefix_size, num_threads);
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < temp_out.size(); i++) {
        temp_out[i] = temp_out[i]/(float)(fft_size);
    }


    //Repeating OFDM syms for each tx antenna with pre and post offset between repetitions of syms
    if (num_reps < num_threads) {
        for (int reps = 0; reps < num_reps; reps++) {
            #pragma omp parallel for num_threads(num_threads)
            for (int tx = 0; tx < num_tx_ants; tx++) {
                //std::cout << "Num threads: " << omp_get_num_threads() << "\n";
                if (reps == 0) {
                    //memcpy((void *)&out_vec[tx][pre_offset + tx*(fft_size + prefix_size)], (const void *)&temp_out[0], fft_size + prefix_size);
                    for (int i = 0; i < fft_size + prefix_size; i++) {
                        out_vec[tx][pre_offset + tx*(fft_size + prefix_size) + i] = temp_out[i];
                    }
                } else {
                    //memcpy((void *)&out_vec[tx][reps*(pre_offset + num_tx_ants*(fft_size + prefix_size) + post_offset) + tx*(fft_size + prefix_size)], (const void *)&temp_out[0], fft_size + prefix_size);
                    for (int i = 0; i < fft_size + prefix_size; i++) {
                        out_vec[tx][reps*(pre_offset + num_tx_ants*(fft_size + prefix_size) + post_offset) + tx*(fft_size + prefix_size) + i] = temp_out[i];
                    }
                }
            }
            #pragma omp barrier
        }
    } else {
        for (int tx = 0; tx < num_tx_ants; tx++) {
            #pragma omp parallel for num_threads(num_threads)
            for (int reps = 0; reps < num_reps; reps++) {
                //std::cout << "Num threads: " << omp_get_num_threads() << "\n";
                if (reps == 0) {
                    //memcpy((void *)&out_vec[tx][pre_offset + tx*(fft_size + prefix_size)], (const void *)&temp_out[0], fft_size + prefix_size);
                    for (int i = 0; i < fft_size + prefix_size; i++) {
                        out_vec[tx][pre_offset + tx*(fft_size + prefix_size) + i] = temp_out[i];
                    }
                } else {
                    //memcpy((void *)&out_vec[tx][reps*(pre_offset + num_tx_ants*(fft_size + prefix_size) + post_offset) + tx*(fft_size + prefix_size)], (const void *)&temp_out[0], fft_size + prefix_size);
                    for (int i = 0; i < fft_size + prefix_size; i++) {
                        out_vec[tx][reps*(pre_offset + num_tx_ants*(fft_size + prefix_size) + post_offset) + tx*(fft_size + prefix_size) + i] = temp_out[i];
                    }
                }
            }
            #pragma omp barrier
        }
    }
    //std::cout << "Sounding frame created...for "<< num_tx_ants << "\n";
    //fftwf_cleanup();

}

void demod_ofdm_symbol(std::complex<float> *in_vec, std::complex<float> *out_vec, int fft_size, int prefix_size) {
    std::complex<float> *temp;
    temp = (std::complex<float> *)malloc((size_t)(fft_size * sizeof(std::complex<float>)));

    //Performing FFT and removing cyclic prefix
    //win_fft(&in_vec[prefix_size], &temp[0], fft_size, -1);
    single_thread_fft(&in_vec[prefix_size], &temp[0], fft_size);

    //Swapping halves
    swap_halves(&temp[0], fft_size);

    //Copying all sub-carrier values except for DC sub-carrier
    memcpy((void *)&out_vec[0], (const void *)&temp[0], ((int)((float)fft_size/(float)2))*sizeof(std::complex<float>));
    //out[fft_size/2] = 0;
    memcpy((void *)&out_vec[(int)((float)fft_size/(float)2)], (const void *)&temp[(int)((float)fft_size/(float)2) + 1], ((int)((float)fft_size/(float)2)-1)*sizeof(std::complex<float>));

    free(temp);
}

void sound_frame(std::vector<std::complex<float>> &pilot_vec, std::vector<std::vector<std::complex<float>>> &in_vec, std::vector<std::vector<std::complex<float>>> &out_vec, int fft_size, int prefix_size, int num_syms, int num_rx_ants, int num_tx_ants, int num_threads) {
    //int num_syms = in_vec[0].size();
    //std::vector<int> fft_size_(num_threads, fft_size), prefix_size_(num_threads, prefix_size);
    out_vec.resize(in_vec.size());
    for (int i = 0; i < in_vec.size(); i++) {
        out_vec[i].resize(num_syms*(fft_size - 1));
    }

    //std::vector<std::thread> thread_vec(num_threads);
    //std::vector<std::complex<float>> *in_ptr(num_threads), *out_ptr(num_threads);
    //omp_set_dynamic(0);
    //omp_set_num_threads(num_threads);


    if (num_syms >= num_threads) {
        //std::cout << "Num syms greater...\n";
        for (int rx = 0; rx < num_rx_ants; rx++) {
            #pragma omp parallel for num_threads(num_threads)
            for (int sym = 0; sym < num_syms; sym += 1) {
                //std::cout << "Sym val: " << sym << "\n";
                //std::cout << "Num threads: " << num_threads << "\n";
                demod_ofdm_symbol(&in_vec[rx][sym*(fft_size + prefix_size)], &out_vec[rx][sym*(fft_size - 1)], fft_size, prefix_size);
                divide_by_vec(&out_vec[rx][(sym)*(fft_size - 1)], &pilot_vec[0], &out_vec[rx][(sym)*(fft_size - 1)], fft_size - 1);
            }
        }
    } else {
        int num_rx_threads = std::max(1,(int)floor(((float)num_rx_ants/(float)(num_syms + num_rx_ants))*num_threads));
        int num_sym_threads = std::max(1,(int)floor((float)num_threads/(float)num_rx_threads));
        int extra_procs = num_threads - (num_rx_threads*num_sym_threads);
        //std::cout << "Num ants greater...\n";
        #pragma omp parallel for num_threads(num_rx_threads)
        for (int rx = 0; rx < num_rx_ants; rx += 1) {
            int offset = 0;
            if (rx < extra_procs) {
                offset = 1;
            }
            #pragma omp parallel for num_threads(num_sym_threads + offset)
            for (int sym = 0; sym < num_syms; sym++) {
                demod_ofdm_symbol(&in_vec[rx][sym*(fft_size + prefix_size)], &out_vec[rx][sym*(fft_size - 1)], fft_size, prefix_size);
                divide_by_vec(&out_vec[rx][(sym)*(fft_size - 1)], &pilot_vec[0], &out_vec[rx][(sym)*(fft_size - 1)], fft_size - 1);
            }
        }
        #pragma omp barrier
    }
    //Averaging multiple vectors for each tx antenna
    if ((int)floor((float)num_syms/(float)num_tx_ants) == 1) {
        return;
    }
    int num_rx_threads = std::max(1,(int)floor(((float)num_rx_ants/(float)(num_tx_ants + num_rx_ants))*num_threads));
    int num_tx_threads = std::max(1,(int)floor((float)num_threads/(float)num_rx_threads));
    int extra_procs = num_threads - (num_rx_threads*num_tx_threads);

    #pragma omp parallel for num_threads(num_rx_threads)
    for (int rx = 0; rx < num_rx_ants; rx += 1) {
        int offset = 0;
        if (rx < extra_procs) {
            offset = 1;
        }
        #pragma omp parallel for num_threads(num_tx_threads + offset)
        for (int tx = 0; tx < num_tx_ants; tx += 1) {
            std::vector<std::complex<float>> temp((fft_size - 1) * (int)floor((float)num_syms/(float)num_tx_ants));
            for (int i = 0; i < num_syms; i += num_tx_ants) {
                memcpy((void *)&temp[(int)floor((float)i/(float)num_tx_ants)*(fft_size - 1)], (const void *)&out_vec[rx][(i + tx)*(fft_size - 1)], (fft_size - 1)*sizeof(std::complex<float>));
            }
            vector_averaging(&temp[0], &temp[0], (int)floor((float)num_syms/(float)num_tx_ants), fft_size - 1);
            memcpy((void *)&out_vec[rx][(tx)*(fft_size - 1)], (const void *)&temp[0], (fft_size - 1)*sizeof(std::complex<float>));
        }
        out_vec[rx].resize(num_tx_ants*(fft_size - 1));
    }
}

/*********************************************************
 * PN Sequence Channel Sounding
 * ******************************************************/


void create_pn_seq_frame(std::vector<std::vector<int>> polys, std::vector<std::vector<std::complex<float>>> &out_vec, int num_tx_ants, int pre_ants, int total_tx_ants, float samp_rate, float max_frame_time, int num_threads) {

    std::vector<std::vector<int>> pn_seq(num_tx_ants);
    int pn_len = (int)pow(2, polys[0].size()-1) - 1;
    std::vector<int> temp_pn;
    if (polys.size() >= total_tx_ants) {
        for (int i = pre_ants; i < pre_ants + num_tx_ants; i++) {
            pn_seq[i - pre_ants] = pn_seq_gen(polys[i], pn_len);
        }
    } else {
        for (int i = pre_ants; i < pre_ants + num_tx_ants; i++) {
            temp_pn = pn_seq_gen(polys[i % polys.size()], pn_len);
            pn_seq[i - pre_ants].resize(pn_len);
            int circshift = i * (int)std::floor((float)pn_len/(float)std::floor((float)total_tx_ants/(float)polys.size()));
            for (int j = 0; j < pn_len; j++) {
                pn_seq[i - pre_ants][j] = temp[(j + circshift) % pn_len];
            }
        }
    }

    int total_frame_samps = (int)std::floor(max_frame_time * samp_rate);
    int num_reps = (int)pow(2, (int)std::floor(log2((float)total_frame_samps/(float)pn_len)));

    wh_matrix wh_mat;
    wh_mat.create_walsh_mat((int)log2(num_reps));

    if (out_vec.size() < num_tx_ants) {
        out_vec.resize(num_tx_ants);
    }
    for (int i = 0; i < num_tx_ants; i++) {
        out_vec[i].resize(pn_len);
    }


    //Repeating PN sequence syms for each tx antenna with pre and post offset between repetitions of syms
    if (num_reps < num_threads) {
        for (int reps = 0; reps < num_reps; reps++) {
            #pragma omp parallel for num_threads(num_threads)
            for (int tx = 0; tx < num_tx_ants; tx++) {
                for (int i = 0; i < pn_len; i++) {
                    out_vec[tx][reps*pn_len + i] = wh_mat.wh_mat[tx][reps] * std::complex<float>(pn_seq[i],0);
                }
            }
            #pragma omp barrier
        }
    } else {
        for (int tx = 0; tx < num_tx_ants; tx++) {
            #pragma omp parallel for num_threads(num_threads)
            for (int reps = 0; reps < num_reps; reps++) {
                //memcpy((void *)&out_vec[tx][reps*(pre_offset + num_tx_ants*(fft_size + prefix_size) + post_offset) + tx*(fft_size + prefix_size)], (const void *)&temp_out[0], fft_size + prefix_size);
                for (int i = 0; i < pn_len; i++) {
                    out_vec[tx][reps*pn_len + i] = wh_mat.wh_mat[tx][reps] * std::complex<float>(pn_seq[i],0);
                }
            }
            #pragma omp barrier
        }
    }

}

void sound_pn_frame(std::vector<std::vector<int>> polys, std::vector<std::vector<std::complex<float>>> &in_vec, std::vector<std::vector<std::complex<float>>> &out_vec, int total_tx_ants, int num_rx_ants, float samp_rate, float max_frame_time, int num_threads) {
    
    int pn_len = (int)pow(2, polys[0].size()-1) - 1;
    std::vector<std::vector<std::complex<float>>> pn_seq(total_tx_ants, std::vector<std::complex<float>>(pn_len));
    std::vector<int> temp_pn;
    if (polys.size() >= total_tx_ants) {
        for (int i = 0; i < total_tx_ants; i++) {
            temp_pn[i] = pn_seq_gen(polys[i], pn_len);
            for (int j = 0; i < pn_len; j++) {
                pn_seq[i][j] = std::complex<float>(temp_pn[i][j], 0);
            }
        }
    } else {
        for (int i = 0; i < total_tx_ants; i++) {
            temp_pn = pn_seq_gen(polys[i % polys.size()], pn_len);
            pn_seq[i].resize(pn_len);
            int circshift = i * (int)std::floor((float)pn_len/(float)std::floor((float)total_tx_ants/(float)polys.size()));
            for (int j = 0; j < pn_len; j++) {
                pn_seq[i][j] = std::complex<float>(temp[(j + circshift) % pn_len], 0);
            }
        }
    }

    int total_frame_samps = (int)std::floor(max_frame_time * samp_rate);
    int num_reps = (int)pow(2, (int)std::floor(log2((float)total_frame_samps/(float)pn_len)));

    wh_matrix wh_mat;
    wh_mat.create_walsh_mat((int)log2(num_reps));
    
    out_vec.resize(in_vec.size());
    for (int i = 0; i < in_vec.size(); i++) {
        out_vec[i].resize(num_reps*total_tx_ants*(2*pn_len - 1));
    }

    //std::vector<std::thread> thread_vec(num_threads);
    //std::vector<std::complex<float>> *in_ptr(num_threads), *out_ptr(num_threads);
    //omp_set_dynamic(0);
    //omp_set_num_threads(num_threads);

    std::vector<std::complex<float>> pn_with_whmat(pn_len);
    if (num_reps >= num_threads) {
        //std::cout << "Num syms greater...\n";
        for (int rx = 0; rx < num_rx_ants; rx++) {
            #pragma omp parallel for num_threads(num_threads)
            for (int sym = 0; sym < num_reps; sym += 1) {
                for (int tx = 0; tx < total_tx_ants; tx++) {
                    for (int j = 0; j < pn_len; j++) {
                        pn_with_whmat[j] = wh_mat.wh_mat[tx][sym] * pn_seq[tx][j];
                    }
                    lin_corr_fft(&in_vec[rx][sym*pn_len], &pn_with_whmat[0], &out_vec[rx][tx*sym*(2*pn_len -1 1)], pn_len, pn_len);
                }
            }
        }
    } else {
        int num_rx_threads = std::max(1,(int)floor(((float)num_rx_ants/(float)(num_reps + num_rx_ants))*num_threads));
        int num_sym_threads = std::max(1,(int)floor((float)num_threads/(float)num_rx_threads));
        int extra_procs = num_threads - (num_rx_threads*num_sym_threads);
        //std::cout << "Num ants greater...\n";
        #pragma omp parallel for num_threads(num_rx_threads)
        for (int rx = 0; rx < num_rx_ants; rx += 1) {
            int offset = 0;
            if (rx < extra_procs) {
                offset = 1;
            }
            #pragma omp parallel for num_threads(num_sym_threads + offset)
            for (int sym = 0; sym < num_reps; sym += 1) {
                for (int tx = 0; tx < total_tx_ants; tx++) {
                    for (int j = 0; j < pn_len; j++) {
                        pn_with_whmat[j] = wh_mat.wh_mat[tx][sym] * pn_seq[tx][j];
                    }
                    lin_corr_fft(&in_vec[rx][sym*pn_len], &pn_with_whmat[0], &out_vec[rx][tx*sym*(2*pn_len - 1)], pn_len, pn_len);
                }
            }
        }
        #pragma omp barrier
    }
    //Averaging multiple vectors for each tx antenna
    if ((int)floor((float)num_syms/(float)total_tx_ants) == 1) {
        return;
    }

    pn_len = 2*pn_len-1;
    int num_rx_threads = std::max(1,(int)floor(((float)num_rx_ants/(float)(total_tx_ants + num_rx_ants))*num_threads));
    int num_tx_threads = std::max(1,(int)floor((float)num_threads/(float)num_rx_threads));
    int extra_procs = num_threads - (num_rx_threads*num_tx_threads);

    #pragma omp parallel for num_threads(num_rx_threads)
    for (int rx = 0; rx < num_rx_ants; rx += 1) {
        int offset = 0;
        if (rx < extra_procs) {
            offset = 1;
        }
        #pragma omp parallel for num_threads(num_tx_threads + offset)
        for (int tx = 0; tx < total_tx_ants; tx += 1) {
            std::vector<std::complex<float>> temp(pn_len * (int)floor((float)num_syms/(float)total_tx_ants));
            for (int i = 0; i < num_syms; i += total_tx_ants) {
                memcpy((void *)&temp[(int)floor((float)i/(float)total_tx_ants)*pn_len], (const void *)&out_vec[rx][(i + tx)*pn_len], pn_len*sizeof(std::complex<float>));
            }
            vector_averaging(&temp[0], &temp[0], (int)floor((float)num_syms/(float)total_tx_ants), pn_len);
            memcpy((void *)&out_vec[rx][(tx)*pn_len], (const void *)&temp[0], pn_len*sizeof(std::complex<float>));
        }
        out_vec[rx].resize(total_tx_ants*pn_len);
    }
}

