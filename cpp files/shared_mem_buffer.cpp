#include "shared_mem_buffer.h"

shared_mem_buffer::shared_mem_buffer() {}

//Creating shared memory with master/slave mode
shared_mem_buffer::shared_mem_buffer(std::string mem, bool ismaster, int size_per_copy_ = 0, int num_ants_ = 0) {
    shared_mem_ptr = new CSharedMemSimple(mem, sizeof(struct frame_buffer));
    buffer = (struct frame_buffer *)shared_mem_ptr->ptr();
    if (ismaster == true) {
        std::cout << "Num workers: " << buffer->num_workers << "\n";
        buffer->size_of_buff = (FFTsize + PREFIXsize)*TXants*RXants*NUMavg + OTHERsize;
        buffer->size_per_copy = size_per_copy_;
        buffer->num_ants = num_ants_;
        if (size_per_copy_ > buffer->size_of_buff) {
            std::cout << "Data cannot fit in memory. Increase shared memory size...\n";
        }
        buffer->rd_offset = -1;
        buffer->wrt_offset = -1;
        is_master = ismaster;
        shared_mem_ptr->set_master_mode();
    } else {
        std::cout << "Num workers increased...\n";
        buffer->num_workers += 1;
    }
}

//Destroying shared memory only if in master mode
shared_mem_buffer::~shared_mem_buffer() {
    if (is_master == true) {
        /*
        while (buffer->num_workers > 0) {
            std::cout << "Num workers left: " << buffer->num_workers << "\n";
            //for (int i = 0; i < 10; i++) {
            //    int j = i;
           // }
        }
        */
        std::cout << "Deleting shared memory...\n";
        delete shared_mem_ptr;
    } else {
        buffer->num_workers -= 1;
//        buffer->size_of_buff = -1;
        std::cout << "Num workers decreased to " << buffer->num_workers << "\n";
    }
}

//Writing to shared memory
void shared_mem_buffer::write_data(std::complex<float> *in) {

    if (buffer->wrt_offset == -1) {
        std::cout << "Writing in shared memory...\n";
        buffer->wrt_offset += 1;
    }

    while (buffer->wrt_offset == buffer->rd_offset) {
        for (int i = 0; i < 10; i++) {
            int j = i;
        }
        std::cout << "Writer Waiting...\n";
    }

    memcpy(&buffer->buffer[buffer->wrt_offset*buffer->size_per_copy], in, buffer->size_per_copy*sizeof(std::complex<float>));
    buffer->wrt_offset += 1;
    if (buffer->wrt_offset >= buffer->num_ants) {
        buffer->wrt_offset = 0;
    }
//    std::cout << "Writer offset iter: " << buffer->wrt_offset << "\n";
}

//Reading from shared memory
void shared_mem_buffer::read_data(std::complex<float> *out) {


    while (buffer->rd_offset + 1 == buffer->wrt_offset || buffer->wrt_offset == -1) {
        //for (int i = 0; i < 10; i++) {
        //    int j = i;
        //}
        //std::cout << "Write pointer: " << buffer->wrt_offset << "\n";
        //std::cout << "Num ants: " << buffer->num_ants << "\n";
        //std::cout << "Reader Waiting...\n";
    }

    //std::cout << "Reader reading...\n";
    memcpy(out, &buffer->buffer[(buffer->rd_offset+1)*buffer->size_per_copy], buffer->size_per_copy*sizeof(std::complex<float>));
    buffer->rd_offset += 1;
    if ((buffer->rd_offset+1) >= buffer->num_ants) {
        buffer->rd_offset = -1;
    }
}

//Returns the number of antennas for which frames are being shared in the shared memory
int shared_mem_buffer::get_num_ants() {
    return buffer->num_ants;
}