#include "chan_sounder_ofdm_mpi.h"

chan_sounder_ofdm_mpi::chan_sounder_ofdm_mpi() {}

chan_sounder_ofdm_mpi::chan_sounder_ofdm_mpi(int fft_size_, int prefix_size_, int num_tx_ants_, int num_rx_ants_, int crank_, int csize_, MPI_Comm *ofdm_comm_) {
    fft_size = fft_size_;
    prefix_size = prefix_size_;
    num_tx_ants = num_tx_ants_;
    num_rx_ants = num_rx_ants_;
    set_fft_size(fft_size_);
    set_cp_size(prefix_size_);
    set_num_tx_ants(num_tx_ants_);
    set_num_rx_ants(num_rx_ants_);
    crank = crank_;
    csize = csize_;
    ofdm_comm = &ofdm_comm_;
}

chan_sounder_ofdm_mpi::~chan_sounder_ofdm_mpi() {}

//Wrapper for function which creates OFDM pilot frame
void chan_sounder_ofdm_mpi::create_pilot_frame(std::vector<std::complex<float>> &pilot_vec, std::vector<std::vector<std::complex<float>>> &out_vec, int pre_offset, int post_offset, int num_reps, int num_local_threads) {
    create_ofdm_sounding_frame(pilot_vec, out_vec, pre_offset, post_offset, num_reps, num_local_threads);
}

//Demodulates OFDM frame for rx antenna data present in the same server
void chan_sounder_ofdm_mpi::demod_ofdm_frame(std::vector<std::complex<float>> &in_vec, std::vector<std::vector<std::complex<float>>> &out_vec) {
    
}
