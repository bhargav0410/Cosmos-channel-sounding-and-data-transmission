#ifndef SHARED_MEM_BUFFER_H
#define SHARED_MEM_BUFFER_H

#include "CSharedMemSimple.hpp"
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <complex.h>
#include <iostream>

#define FFTsize 1024
#define PREFIXsize 64
#define TXants 16
#define RXants 16
#define NUMavg 100
#define OTHERsize 0

struct frame_buffer {
    std::complex<float> buffer[(FFTsize + PREFIXsize)*TXants*RXants*NUMavg + OTHERsize];
    int size_per_copy, rd_offset = -1, wrt_offset = -1, size_of_buff, num_ants, num_workers = 0;
};

class shared_mem_buffer {

public:
    shared_mem_buffer();
    shared_mem_buffer(std::string, bool, int, int);
    ~shared_mem_buffer();
    void write_data(std::complex<float> *);
    void read_data(std::complex<float> *);
    int get_num_ants();

private:
    frame_buffer *buffer;
    CSharedMemSimple* shared_mem_ptr;
    bool is_master = false;
    int actual_num_syms;

};

#endif