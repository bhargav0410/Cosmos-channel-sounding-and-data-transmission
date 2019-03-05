#include "chan_sound_pn_seq_gpu.cuh"
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <random>
#include <iostream>

int main(int argc, char *argv[]) {
    std::string file;
    int pn_length, num_tx_seqs, num_rx_ants;

    //Taking input arguments for PN sequence file and length and number of sequences
    if (argc > 1) {
        file = argv[1];
    }
    if (argc > 2) {
        pn_length = atoi(argv[2]);
    }
    if (argc > 3) {
        num_tx_seqs = atoi(argv[3]);
    }
    if (argc > 4) {
        num_rx_ants = atoi(argv[4]);
    }

    //Taking input sequences from file
    std::ifstream infile(file.c_str(), std::ifstream::binary);
    std::complex<float> *temp_arr;
    temp_arr = (std::complex<float> *)malloc(pn_length*num_tx_seqs*sizeof(*temp_arr));
    infile.read((char *)temp_arr, pn_length*num_tx_seqs*sizeof(*temp_arr));

    //Adding noise to input PN sequences
    std::default_random_engine generator;
    std::normal_distribution<float> distribution(0.0,0.5);

    std::vector<std::complex<float>> noise(pn_length);
    for (int i = 0; i < pn_length; i++) {
        noise[i] = std::complex<float>(distribution(generator), distribution(generator));
    }

    std::vector<std::vector<std::complex<float>>> rec_seq(num_rx_ants, std::vector<std::complex<float>>(pn_length));
    //Adding all the PN sequences together and then adding noise
    for (int tx = 0; tx < num_tx_seqs; tx++) {
        for (int rx = 0; rx < num_rx_ants; rx++) {
            for (int i = 0; i < pn_length; i++) {
                rec_seq[rx][i] += temp_arr[tx*pn_length + i];
            }
        }
    }
    for (int rx = 0; rx < num_rx_ants; rx++) {
        for (int i = 0; i < pn_length; i++) {
            rec_seq[rx][i] = rec_seq[rx][i]/(std::complex<float>((float)num_tx_seqs, 0.0));
            rec_seq[rx][i] += noise[i];
        }
    }

    //Printing the received and noisy signal
    for (int rx = 0; rx < num_rx_ants; rx++) {
        std::cout << "RX " << rx + 1 << "\n";
        for (int i = 0; i < pn_length; i++) {
            std::cout << rec_seq[rx][i] << " ";
        }
        std::cout << "\n";
    }

    //Creating object for GPU correlation
    chan_sound_pn_seq_gpu pn_seq_gpu(255, 1, 1, 0);

    pn_seq_gpu.get_pn_seq_from_file(file, num_tx_seqs, pn_length);

    //Output array
    std::vector<std::vector<std::vector<std::complex<float>>>> out_array;

    //Correlating with PN sequence
    pn_seq_gpu.correlate_with_pn_seq(rec_seq, out_array);

    //Showing output array
    for (int tx = 0; tx < num_tx_seqs; tx++) {
        std::cout << "TX " << tx + 1 << "\n";
        for (int rx = 0; rx < num_rx_ants; rx++) {
            std::cout << "RX " << rx + 1 << "\n";
            for (int i = 0; i < out_array[tx][rx].size(); i++) {
                std::cout << out_array[tx][rx][i] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n-------------------------\n\n";
    }

    free(temp_arr);
    return 0;
}