#ifndef CHAN_SOUNDER_OFDM_MPI_H
#define CHAN_SOUNDER_OFDM_MPI_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <vector>
#include <complex>
#include <cmath>
#include <thread>
#include <omp.h>
#include <mpi.h>
#include "chan_sounder.h"

class chan_sounder_ofdm_mpi : public chan_sounder {

public:
    chan_sounder_ofdm_mpi();
    chan_sounder_ofdm_mpi(int, int, int, int, int, int, MPI_Comm *);
    ~chan_sounder_ofdm_mpi();
    //void load_balancing_mpi(int *, int *, int, int);
    void create_pilot_frame(std::vector<std::complex<float>> &, std::vector<std::vector<std::complex<float>>> &, int, int, int, int);

private:
    fft_size, prefix_size, num_tx_ants, num_rx_ants, crank, csize;
    MPI_comm ofdm_comm;

}

#endif