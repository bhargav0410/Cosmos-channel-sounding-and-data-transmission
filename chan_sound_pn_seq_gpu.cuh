#ifndef CHAN_SOUND_PN_SEQ_GPU_CUH
#define CHAN_SOUND_PN_SEQ_GPU_CUH

#include <string>
#include <cstdlib>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <fstream>

class chan_sound_pn_seq_gpu {

public:
    chan_sound_pn_seq_gpu() {}

    chan_sound_pn_seq_gpu(int _pn_len, int _num_tx_seqs, int _num_rx_ants, int _thres) {
        pn_length = _pn_len;
        num_tx_seqs = _num_tx_seqs;
        num_rx_ants = _num_rx_ants;
        thres = _thres;
    }

    ~chan_sound_pn_seq_gpu() {}

    void get_pn_seq_from_file(std::string file, int seqs_in_file) {

    }

    void 

private:
    int pn_length, num_tx_seqs, num_rx_ants, thres;

};

#endif