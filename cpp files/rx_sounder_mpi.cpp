/*
Copyright (c) 2018, WINLAB, Rutgers University, USA
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <uhd/utils/thread.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <uhd/types/time_spec.hpp>
#include <uhd/transport/udp_simple.hpp>
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
#include <sys/socket.h>
#include <string>
#include <omp.h>
#include <mpi.h>
#include "chan_sounder.h"
#include "shared_mem_buffer.h"

//#include <ofdm_mpi.h>

#define pi 3.141592654
#define c 299792458.0

namespace po = boost::program_options;
using boost::asio::ip::udp;
/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}
size_t len = 0;
void handler(
  const boost::system::error_code& error, // Result of operation.
  size_t length          // Number of bytes received.
  ) {
	  len = length;
  } 

  
/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char *argv[]){
    //uhd::set_thread_priority_safe();
	//Initialzing MPI
    int provided;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
	//getting size and rank
	int gsize, grank;
	MPI_Comm_size(MPI_COMM_WORLD, &gsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    //variables to be set by po
    std::string dist_args, args, arg_ip_part, ant, subdev, ref, pps, otw, cpu, channel_list, file, addr, remote_addr;
    double rate, freq, gain, bw;
    float delay, thres, sounding_time;
	int port, num_times, chan_per_radio;
	size_t nsamps; 
	bool pn_seq, same;
	int start_num, num_usrp;
	int fft_size, prefix_size, num_reps, num_threads, num_tx_ants;
    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
		//("start_of_addr", po::value<std::string>(&arg_ip_part)->default_value(""), "common part of IP address")
		//("start_num_of_addr", po::value<int>(&start_num)->default_value(1), "start number of IP address")
		//("nusrp", po::value<int>(&num_usrp)->default_value(1), "Number of USRPs to be used")
        ("args", po::value<std::string>(&args)->default_value(""), "IP address of USRP devices")
        ("rate", po::value<double>(&rate), "rate of outgoing samples")
        ("freq", po::value<double>(&freq), "RF center frequency in Hz")
        ("gain", po::value<double>(&gain), "gain for the RF chain")
		("ant", po::value<std::string>(&ant), "antenna selection")
        ("subdev", po::value<std::string>(&subdev), "subdevice specification")
        ("bw", po::value<double>(&bw), "analog frontend filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo, gpsdo)")
        ("sync", po::value<std::string>(&pps), "Time sync source (now, pps, mimo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode (sc16 or sc8)")
		("cpu", po::value<std::string>(&cpu)->default_value("fc32"), "specify the cpu sample mode (fc32 or sc16)")
        ("chan_per_usrp", po::value<int>(&chan_per_radio)->default_value(1), "number of antennas to be used per USRP")
		//("channels", po::value<std::string>(&channel_list)->default_value("0"), "which channels to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("int-n", "tune USRP with integer-N tuning")
		("delay", po::value<float>(&delay)->default_value(0), "Delay value in stream time spec")
		("sounder_resolution", po::value<int>(&fft_size)->default_value(64), "Resolution of channel sounder (in powers of 2)")
		("prefix_size", po::value<int>(&prefix_size)->default_value(16), "Size of cyclic prefix (greater than expected channel spread and less than sounder resolution size)")
		("num_tx_ants", po::value<int>(&num_tx_ants)->default_value(1), "Number of transmitting antennas for which channel sounding is being performed")
        ("avg_window", po::value<int>(&num_reps)->default_value(1), "Number of sounding symbols per Tx antenna to be averaged (minimum is 1 which means no averaging)")
		("num_threads", po::value<int>(&num_threads)->default_value(1), "Number of threads per CPU for packet finding")
        ("threshold", po::value<float>(&thres)->default_value(0.1), "Threshold for finding packet")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")){
        std::cout << boost::format("UHD RX to shared memory %s") % desc << std::endl;
        return ~0;
	}

	std::vector<std::string> arg_strings;
    boost::split(arg_strings, args, boost::is_any_of(","));
	num_usrp = arg_strings.size();

    //create a usrp device

	std::cout << "Number of USRPs selected = " << num_usrp << std::endl;
	std::vector<int> size_of_proc_data(gsize), displ(gsize);
	if (num_usrp > 1) {
		
		int displ_of_proc = 0;
		for (int i = 0; i < gsize; i++) {
			displ[i] = displ_of_proc;
			size_of_proc_data[i] = (int)floor((float)num_usrp/(float)gsize);
			if (i < (int)num_usrp % gsize) {
				size_of_proc_data[i] += 1;
			}
			displ_of_proc += size_of_proc_data[i];
		}
		std::cout << "Start point: " << displ[grank] << std::endl;
		std::cout << "Length: " << size_of_proc_data[grank] << std::endl;
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i = displ[grank]; i < displ[grank] + size_of_proc_data[grank]; i++) {
			dist_args += "addr" + std::to_string(i - displ[grank]) + "=" + arg_strings[i] + ",";
			//std::cout << args;
		}
	} else {
		size_of_proc_data[0] = 1;
		displ[0] = 0;
		dist_args = dist_args + "addr0=" + arg_strings[0];
	}

	std::cout << "Address used : " << dist_args;

    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % dist_args << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(dist_args);

	std::vector<size_t> channel_nums;
	for (int i = 0; i < size_of_proc_data[grank]; i++) {
		for (int ch = 0; ch < chan_per_radio; ch++) {
			channel_nums.push_back(ch + i*chan_per_radio);
		}
	}
    std::cout << "Total no. of channels: " << channel_nums.size() << "\n";


    //Lock mboard clocks
    usrp->set_clock_source(ref);

     //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev")) usrp->set_rx_subdev_spec(subdev);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    //set the sample rate
    if (not vm.count("rate")){
        std::cerr << "Please specify the sample rate with --rate" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rate/1e6) << std::endl;
    usrp->set_rx_rate(rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...") % (usrp->get_rx_rate()/1e6) << std::endl << std::endl;

	for (int ch = 0; ch < channel_nums.size(); ch++) {
		//set the center frequency
		if (not vm.count("freq")){
			std::cerr << "Please specify the center frequency with --freq" << std::endl;
			return EXIT_FAILURE;
		}
		
		uhd::tune_request_t tune_request(freq, ch);
		if(vm.count("int-n")) tune_request.args = uhd::device_addr_t("mode_n=integer");
		usrp->set_rx_freq(tune_request, ch);
		std::cout << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq()/1e6) << std::endl << std::endl;

		//set the rf gain
		if (vm.count("gain")){
			std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
			usrp->set_rx_gain(gain, ch);
			std::cout << boost::format("Actual RX Gain: %f dB...") % usrp->get_rx_gain() << std::endl << std::endl;
		}

		//set the analog frontend filter bandwidth
		if (vm.count("bw")){
			std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (bw/1e6) << std::endl;
			usrp->set_rx_bandwidth(bw, ch);
			std::cout << boost::format("Actual RX Bandwidth: %f MHz...") % (usrp->get_rx_bandwidth()/1e6) << std::endl << std::endl;
		}

		//set the antenna
		if (vm.count("ant")) {
			std::cout << boost::format("Setting RX Antenna: %f...") % (ant) << std::endl;
			usrp->set_rx_antenna(ant, ch);
			std::cout << boost::format("Actual RX Antenna: %f...") % (usrp->get_rx_antenna(ch)) << std::endl << std::endl;
		}
	}

    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
    if (pps == "now"){
        //This is not a true time lock, the devices will be off by a few RTT.
        //Rather, this is just to allow for demonstration of the code below.
        usrp->set_time_now(uhd::time_spec_t(0.0));
    }
    else if (pps == "pps"){
        usrp->set_time_source("external");
        usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
        boost::this_thread::sleep(boost::posix_time::seconds(1)); //wait for pps sync pulse
    }
    else if (pps == "mimo"){
        UHD_ASSERT_THROW(usrp->get_num_mboards() == 2);

        //make mboard 1 a slave over the MIMO Cable
        usrp->set_clock_source("mimo", 1);
        usrp->set_time_source("mimo", 1);

        //set time on the master (mboard 0)
        usrp->set_time_now(uhd::time_spec_t(0.0), 0);

        //sleep a bit while the slave locks its time to the master
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
	
	
    boost::this_thread::sleep(boost::posix_time::seconds(1)); //allow for some setup time

	std::vector<int> pn_buff;
	std::vector<int> poly = {1,0,0,0,1,1,1,0,1};
	pn_buff = pn_seq_gen(poly, 255);


	std::vector<std::complex<float>> pilots(fft_size - 1);
	srand(0);
	//std::cout << "Creating pilot vector...\n";
	for (int i = 0; i < fft_size - 1; i++) {
		pilots[i] = std::complex<float>((float)pn_buff[i % pn_buff.size()], 0);//(float)0.707 * std::complex<float>((float)(2*(rand() % 2) - 1), (float)(2*(rand() % 2) - 1));
	}
	std::vector<std::vector<std::complex<float>>> out_vec;
	

	std::vector<std::vector<std::complex<float> > > buff1(
        channel_nums.size(), std::vector<std::complex<float>>(2*(pn_buff.size() + (fft_size + prefix_size)*num_tx_ants*num_reps))
    );
	std::vector<std::vector<std::complex<float> > > buff2(
        channel_nums.size(), std::vector<std::complex<float>>(pn_buff.size() + (fft_size + prefix_size)*num_tx_ants*num_reps)
    );
    std::vector<std::vector<std::complex<float> > > final_buff1(
        channel_nums.size(), std::vector<std::complex<float>>((fft_size + prefix_size)*num_tx_ants*num_reps)
    );
    std::vector<std::vector<std::complex<float> > > final_buff2(
        channel_nums.size(), std::vector<std::complex<float>>((fft_size + prefix_size)*num_tx_ants*num_reps)
    );
    std::vector<std::vector<std::complex<float> > > final_out;
    std::vector<std::complex<float> > global_out(size_of_proc_data[grank]*chan_per_radio*num_tx_ants*(fft_size - 1));

    shared_mem_buffer shmem("rx_shmem", true, final_buff1[0].size(), channel_nums.size());

	MPI_Barrier(MPI_COMM_WORLD);

    std::vector<std::complex<float> *> buff_ptrs1;
    for (size_t i = 0; i < buff1.size(); i++) buff_ptrs1.push_back(&buff1[i].front());
	std::vector<std::complex<float> *> buff_ptrs2;
    for (size_t i = 0; i < buff2.size(); i++) buff_ptrs2.push_back(&buff2[i].front());

    uhd::stream_args_t stream_args("fc32",otw);
	stream_args.channels = channel_nums;
	uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);
	uhd::rx_metadata_t md;
	
	uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
	stream_cmd.stream_now = false;
	stream_cmd.time_spec = uhd::time_spec_t(2 + delay);
	double timeout = delay + 0.5;
	rx_stream->issue_stream_cmd(stream_cmd);
	
	std::signal(SIGINT, &sig_int_handler);
	std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

	//md.start_of_burst = false;
	//md.end_of_burst   = false;
	//md.has_time_spec = true;
	std::cout << "Starting RX streamer...\n";

    //omp_lock_t lock1, lock2;
    int buff1_flag = 0, buff2_flag = 0, num_rx_samps1 = 0, num_rx_samps2 = 0, buff1_inc = 0, buff2_inc = 0;
    int offset = -1, buff1_found = 0, buff2_found = 0;
 //   int write_unlock1 = -1, read_unlock1 = -1, write_unlock2 = -1, read_unlock2 = -1;
    //omp_init_lock(&lock1);
    //omp_init_lock(&lock2);
    clock_t start, finish, pn_start, pn_finish;

    start = clock();
	while(not stop_signal_called) {
        
		//receive the entire contents to the buffer
        #pragma omp parallel num_threads(2) shared(buff1_flag, buff2_flag)
        {
            #pragma omp single nowait
            {
                #pragma omp task
                {
                    //Waiting for threads to finish channel sounding process
//                    while (buff1_flag == 1) {}
                        num_rx_samps1 = rx_stream->recv(buff_ptrs1, (int)buff1[0].size(), md, timeout);
                        //std::cout << "Num samps recvd: " << num_rx_samps1 << "\n";
                        stream_cmd.stream_now = true;
                        timeout = 0.5;
                        
                        //num_rx_samps2 = rx_stream->recv(buff_ptrs2, (int)buff2[0].size(), md, timeout);
                        
                        //while (buff2_flag == 0) {}

                        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
                            offset = -1;
                        }
                    if (offset == -1) {
                        //omp_set_dynamic(0);
                        //omp_set_num_threads(num_threads - 1);
                        pn_start = clock();
                        offset = find_pn_seq(&buff1[0][0], &pn_buff[0], buff1[0].size(), pn_buff.size(), thres, num_threads);
                        pn_finish = clock();
                        std::cout << "Packet finding time: " << (float)(pn_finish - pn_start)/(float)CLOCKS_PER_SEC;
                        if (offset >= 0) {
                            while (buff1_flag == 1) {
                            for (int i = 0; i < 10; i++) {
                                    int j = i;
                                }
                            }
                                for (int i = 0; i < buff1.size(); i++) {
                                    for (int j = offset + pn_buff.size(); j < offset + pn_buff.size() + final_buff1[i].size(); j++) {
                                        final_buff1[i][j - (offset + pn_buff.size())] = buff1[i][j];
                                    }
                                }
                                buff1_flag = 1;
                            
                            
                            while (buff2_flag == 1) {
                                for (int i = 0; i < 10; i++) {
                                    int j = i;
                                }
                            }
                                int offset_for_buff_2 = offset + 2*pn_buff.size() + final_buff1[0].size();
                                int num_samps_left = (int)buff1[0].size() - offset_for_buff_2;
                                if (num_samps_left > 0) {
                                    for (int i = 0; i < buff1.size(); i++) {
                                        for (int j = offset_for_buff_2; j < buff1[i].size(); j++) {
                                            final_buff2[i][j - offset_for_buff_2] = buff1[i][j];
                                        }
                                    }
                                    int num_samps_to_recv = (int)final_buff2[0].size() - num_samps_left;
                                    num_rx_samps2 = rx_stream->recv(buff_ptrs2, num_samps_to_recv, md, timeout);
                                    for (int i = 0; i < buff2.size(); i++) {
                                        for (int j = 0; j < num_samps_to_recv; j++) {
                                            final_buff2[i][j + num_samps_left] = buff2[i][j];
                                        }
                                    }
                                } else {
                                    offset_for_buff_2 = offset + pn_buff.size() + final_buff1[0].size();
                                    num_samps_left = (int)buff1[0].size() - offset_for_buff_2;
                                    int num_samps_to_recv = pn_buff.size() - num_samps_left + final_buff2[0].size();
                                    num_rx_samps2 = rx_stream->recv(buff_ptrs2, num_samps_to_recv, md, timeout);
                                    for (int i = 0; i < buff2.size(); i++) {
                                        for (int j = 0; j < final_buff2[i].size(); j++) {
                                            final_buff2[i][j] = buff2[i][j + pn_buff.size() - num_samps_left];
                                        }
                                    }
                                }
                                buff2_flag = 1;
                            
                        }
                    } else {
                        while (buff1_flag == 1) {
                            for (int i = 0; i < 10; i++) {
                                int j = i;
                            }
                        }
                            
                            for (int i = 0; i < final_buff1.size(); i++) {
                                for (int j = 0; j < final_buff1[0].size(); j++) {
                                    final_buff1[i][j] = buff1[i][j + pn_buff.size()];
                                }
                            }
                            buff1_flag = 1;
                        

                        while (buff2_flag == 1) {
                            for (int i = 0; i < 10; i++) {
                                int j = i;
                            }
                        }
                        
                            for (int i = 0; i < final_buff2.size(); i++) {
                                for (int j = 0; j < final_buff2[0].size(); j++) {
                                    final_buff2[i][j] = buff1[i][j + 2*pn_buff.size() + final_buff1[0].size()];
                                }
                            }
                            buff2_flag = 1;
                        

                    }
                    //buff1_flag = 0;
                    
                    
                    //Waiting for threads to receive data
                    //while (buff2_flag == 0) {}
                    /*
                    if (offset == -1) {
                        //omp_set_dynamic(0);
                        //omp_set_num_threads(num_threads - 1);
                        offset = find_pn_seq(&buff2[0][0], &pn_buff[0], buff2[0].size(), pn_buff.size(), thres);
                        if (offset == -1) {
                            //std::cout << "Cannot find packet. Please lower threshold...\n";
                        } else {
                            //std::cout << "Offset: " << offset << "\n";
                            for (int i = 0; i < buff2.size(); i++) {
                                for (int j = offset + pn_buff.size(); j < buff2[i].size(); j++) {
                                    final_buff[i][j - (offset + pn_buff.size())] = buff2[i][j];
                                }
                            }
                            buff2_found = 1;
                        }
                    } else {
                        if (buff1_found == 1) {
                            for (int i = 0; i < buff2.size(); i++) {
                                for (int j = 0; j < offset; j++) {
                                    final_buff[i][j + (buff2[i].size() - offset - pn_buff.size())] = buff2[i][j];
                                }
                            }
                        } else {
                            for (int i = 0; i < buff2.size(); i++) {
                                for (int j = offset + pn_buff.size(); j < buff2[i].size(); j++) {
                                    final_buff[i][j - (offset + pn_buff.size())] = buff2[i][j];
                                }
                            }
                        }
                    }
                    */
                   // buff2_flag = 0;
                    //std::cout << "Final Offset: " << offset << "\n";
                }
                #pragma omp task
                {
                    while (buff1_flag == 0) {
                        for (int i = 0; i < 10; i++) {
                            int j = i;
                        }
                    }
                    for (int i = 0; i < final_buff1.size(); i++) {
                        shmem.write_data(&final_buff1[i][0]);
                    }
                    buff1_flag = 0;
                    
                    while (buff2_flag == 0) {
                        for (int i = 0; i < 10; i++) {
                            int j = i;
                        }
                    }
                    for (int i = 0; i < final_buff2.size(); i++) {
                        shmem.write_data(&final_buff2[i][0]);
                    }
                    buff2_flag = 0;
                    

                   // std::cout << "Starting second task...\n";
                //    while (buff2_flag == 1) {}
                    //std::cout << "Got the offset...\n";
                    /*
                    while (buff1_flag == 0);
                        pn_start = clock();

                        //std::cout << "Starting sounding...\n";
                        //Performing channel sounding on final synced buffer
                        sound_frame(pilots, final_buff1, final_out, fft_size, prefix_size, num_tx_ants*num_reps, (int)channel_nums.size(), num_tx_ants, num_threads - 1);
                        buff1_flag = 0;

                        for (int i = displ[grank]*chan_per_radio; i < size_of_proc_data[grank]*chan_per_radio; i++) {
                            for (int j = 0; j < final_out[0].size(); j++) {
                                global_out[i*(final_out[0].size()) + j] = final_out[i][j];
                            }
                        }
                        //std::cout << "Copied part of global matrix. Now copying other parts of matrix using MPI...\n";
                        if (grank == 0) {
                            MPI_Gather(MPI_IN_PLACE, final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, (void *)&global_out[displ[grank]*chan_per_radio*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, 0, MPI_COMM_WORLD);
                        } else {
                            MPI_Gather((void *)&global_out[displ[grank]*chan_per_radio*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, (void *)&global_out[displ[grank]*chan_per_radio*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, 0, MPI_COMM_WORLD);
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        pn_finish = clock();
                        std::cout << "Sounding time: " << (float)(pn_finish - pn_start)/(float)CLOCKS_PER_SEC << "\n";

                    while (buff2_flag == 0);

                        pn_start = clock();
                        sound_frame(pilots, final_buff2, final_out, fft_size, prefix_size, num_tx_ants*num_reps, (int)channel_nums.size(), num_tx_ants, num_threads - 1);
                        buff2_flag = 0;

                        //std::cout << "Frame sounded...\n";
                        //std::cout << "Per rx antenna frame size: " << (float)final_out[0].size()/(float)(fft_size - 1) << "\n";
                        for (int i = displ[grank]*chan_per_radio; i < size_of_proc_data[grank]*chan_per_radio; i++) {
                            for (int j = 0; j < final_out[0].size(); j++) {
                                global_out[i*(final_out[0].size()) + j] = final_out[i][j];
                            }
                        }
                        //std::cout << "Copied part of global matrix. Now copying other parts of matrix using MPI...\n";
                        if (grank == 0) {
                            MPI_Gather(MPI_IN_PLACE, final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, (void *)&global_out[displ[grank]*chan_per_radio*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, 0, MPI_COMM_WORLD);
                        } else {
                            MPI_Gather((void *)&global_out[displ[grank]*chan_per_radio*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, (void *)&global_out[displ[grank]*chan_per_radio*final_out[0].size()], final_out[0].size()*size_of_proc_data[grank]*chan_per_radio, MPI_COMPLEX, 0, MPI_COMM_WORLD);
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        pn_finish = clock();
                        std::cout << "Sounding time: " << (float)(pn_finish - pn_start)/(float)CLOCKS_PER_SEC << "\n";
                    */
                    
                }
            }
        }
        //#pragma omp barrier

        finish = clock();
        //Check for stop signal in case of timed channel sounding
        /*
        if ((int)sounding_time != -1) {
            //std::cout << "Time: " << (float)(finish - start)/(float)CLOCKS_PER_SEC << "\n";
            if (((float)(finish - start))/(float)CLOCKS_PER_SEC >= sounding_time + delay + 2) {
                stop_signal_called = true;
            }
        }
        */

	}
	
    rx_stream->issue_stream_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);

    /*
    std::ofstream outfile;
    std::string outfilename = "PN_Seq_used.dat";
    outfile.open(outfilename.c_str(), std::ofstream::binary);
    outfile.write((const char*)&pn_buff.front(), pn_buff.size()*sizeof(int));
    if (outfile.is_open()) {
        outfile.close();
    }
    outfilename = "global_out.dat";
    outfile.open(outfilename.c_str(), std::ofstream::binary);
    outfile.write((const char*)&global_out.front(), global_out.size()*sizeof(std::complex<float>));
    if (outfile.is_open()) {
        outfile.close();
    }
    */
    std::ofstream outfile;
    for (int i = displ[grank]*chan_per_radio; i < (displ[grank] + size_of_proc_data[grank])*chan_per_radio; i++) {
        std::string outfilename = "ch_" + boost::lexical_cast<std::string>(i) + "_binary";
        outfile.open(outfilename.c_str(), std::ofstream::binary);
        outfile.write((const char*)&final_buff1[i].front(), final_buff1[i].size()*sizeof(std::complex<float>));
        outfile.write((const char*)&final_buff2[i].front(), final_buff2[i].size()*sizeof(std::complex<float>));
        if (outfile.is_open()) {
            outfile.close();
        }
        outfilename = "ch_dumped_" + boost::lexical_cast<std::string>(i) + "_binary";
        std::cout << outfilename << "\n";
        outfile.open(outfilename.c_str(), std::ofstream::binary);
        outfile.write((const char*)&buff1[i].front(), buff1[i].size()*sizeof(std::complex<float>));
        //outfile.write((const char*)&final_buff2[i].front(), final_buff2[i].size()*sizeof(std::complex<float>));
        if (outfile.is_open()) {
            outfile.close();
        }
    }
    
	
	//send a mini EOB packet
	//md.end_of_burst = true;
    //finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

	MPI_Finalize();
    return EXIT_SUCCESS;
}


