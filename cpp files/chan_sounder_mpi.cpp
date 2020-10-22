#include <boost/program_options.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <csignal>
#include <ctime>
#include <fstream>
#include <cmath>
#include <thread>
#include <string>
#include <omp.h>
#include <mpi.h>
#include "chan_sounder.h"
#include "shared_mem_buffer.h"

namespace po = boost::program_options;

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

int main (int argc, char *argv[]) {

    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
	//getting size and rank
	int gsize, grank;
	MPI_Comm_size(MPI_COMM_WORLD, &gsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    //variables to be set by po
    std::string file;
    double rate, freq, gain, bw;
    float delay, thres, sounding_time;
	size_t nsamps; 
	bool pn_seq, same;
	int fft_size, prefix_size, num_reps, num_threads, num_tx_ants, num_rx_ants;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
		("sounder_resolution", po::value<int>(&fft_size)->default_value(64), "Resolution of channel sounder (in powers of 2)")
		("prefix_size", po::value<int>(&prefix_size)->default_value(16), "Size of cyclic prefix (greater than expected channel spread and less than sounder resolution size)")
		("num_tx_ants", po::value<int>(&num_tx_ants)->default_value(1), "Number of transmitting antennas for which channel sounding is being performed")
        //("num_rx_ants", po::value<int>(&num_rx_ants)->default_value(1), "Total number of receiving antennas for which channel sounding is being performed in all CPUs")
        ("avg_window", po::value<int>(&num_reps)->default_value(1), "Number of sounding symbols per Tx antenna to be averaged (minimum is 1 which means no averaging)")
		("num_threads", po::value<int>(&num_threads)->default_value(1), "Number of threads per CPU")
        //("sound_for", po::value<float>(&sounding_time)->default_value(-1), "Time to perform channel sounding for. (Put -1 for indefinite)")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")){
        std::cout << boost::format("Channel sounding %s") % desc << std::endl;
        return ~0;
	}

    shared_mem_buffer shmem("rx_shmem", false, 0, 0);
    num_rx_ants = shmem.get_num_ants();
    std::vector<int> displ(gsize, 0);
    std::vector<int> size_of_proc_data(gsize);
    size_of_proc_data[grank] = num_rx_ants;

    MPI_Barrier(MPI_COMM_WORLD);

    int total_rx_ants = 0;
    for (int i = 0; i < gsize; i++) {
        MPI_Bcast(&size_of_proc_data[i], 1, MPI_INT, i, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        for (int j = i + 1; j < gsize; j++) {
            displ[j] += size_of_proc_data[i];
        }
        total_rx_ants += size_of_proc_data[i];
    }
    

    std::vector<std::vector<std::complex<float> > > final_buff1(
        num_rx_ants, std::vector<std::complex<float>>((fft_size + prefix_size)*num_tx_ants*num_reps)
    );
    std::vector<std::vector<std::complex<float> > > final_out;
    std::vector<std::complex<float> > global_out(total_rx_ants*num_tx_ants*(fft_size - 1));

	std::vector<int> pn_buff;
	std::vector<int> poly = {1,0,0,0,1,1,1,0,1};
	pn_buff = pn_seq_gen(poly, 255);
	std::vector<std::complex<float>> pilots(fft_size - 1);
	srand(0);
	//std::cout << "Creating pilot vector...\n";
	for (int i = 0; i < fft_size - 1; i++) {
		pilots[i] = std::complex<float>((float)pn_buff[i % pn_buff.size()], 0);//(float)0.707 * std::complex<float>((float)(2*(rand() % 2) - 1), (float)(2*(rand() % 2) - 1));
	}
    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop sounding..." << std::endl;

    clock_t start, finish;
    while (not stop_signal_called) {

        start = clock();
        for (int i = 0; i < num_rx_ants; i++) {
            shmem.read_data(&final_buff1[i][0]);
        }

        sound_frame(pilots, final_buff1, final_out, fft_size, prefix_size, num_tx_ants*num_reps, num_rx_ants, num_tx_ants, num_threads);

        for (int i = displ[grank]; i < size_of_proc_data[grank]; i++) {
            for (int j = 0; j < final_out[0].size(); j++) {
                global_out[i*(final_out[0].size()) + j] = final_out[i][j];
            }
        }
        //std::cout << "Copied part of global matrix. Now copying other parts of matrix using MPI...\n";
        if (grank == 0) {
            MPI_Gather(MPI_IN_PLACE, final_out[0].size()*num_rx_ants, MPI_COMPLEX, (void *)&global_out[displ[grank]*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank], MPI_COMPLEX, 0, MPI_COMM_WORLD);
        } else {
            MPI_Gather((void *)&global_out[displ[grank]*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank], MPI_COMPLEX, (void *)&global_out[displ[grank]*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank], MPI_COMPLEX, 0, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        finish = clock();

        if (grank == 0) {
            std::cout << "Processing time: " << (float)(finish - start)/(float)CLOCKS_PER_SEC << "\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    std::ofstream outfile;
    if (grank == 0) {
        std::string outfilename = "global_out.dat";
        outfile.open(outfilename.c_str(), std::ofstream::binary);
        outfile.write((const char*)&global_out.front(), global_out.size()*sizeof(std::complex<float>));
        if (outfile.is_open()) {
            outfile.close();
        }
    }

    std::cout << "Sounding done...\n";
	MPI_Finalize();
    return EXIT_SUCCESS;

}