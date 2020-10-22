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
#include <mpi.h>
#include "chan_sounder.h"

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
    float delay;
	int port, start ,end, num_times, chan_per_radio;
	size_t nsamps; 
	bool pn_seq, same;
	int start_num, num_usrp;
	int fft_size, prefix_size, num_reps, num_threads;
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
        ("pps", po::value<std::string>(&pps), "PPS source (internal, external, mimo, gpsdo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode (sc16 or sc8)")
		("cpu", po::value<std::string>(&cpu)->default_value("fc32"), "specify the cpu sample mode (fc32 or sc16)")
        ("chan_per_usrp", po::value<int>(&chan_per_radio)->default_value(1), "number of antennas to be used per USRP")
		//("channels", po::value<std::string>(&channel_list)->default_value("0"), "which channels to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("int-n", "tune USRP with integer-N tuning")
		("delay", po::value<float>(&delay)->default_value(0), "Delay value in stream time spec")
		("sounder_resolution", po::value<int>(&fft_size)->default_value(64), "Resolution of channel sounder (in powers of 2)")
		("prefix_size", po::value<int>(&prefix_size)->default_value(16), "Size of cyclic prefix (greater than expected channel spread and less than sounder resolution size)")
		("avg_window", po::value<int>(&num_reps)->default_value(1), "Number of sounding symbols per Tx antenna to be averaged (minimum is 1 which means no averaging)")
		("num_threads", po::value<int>(&num_threads)->default_value(1), "Number of threads per CPU")	
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")){
        std::cout << boost::format("UHD TX Waveforms %s") % desc << std::endl;
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
	/*
    //detect which channels to use
    std::vector<std::string> channel_strings;
    std::vector<size_t> channel_nums;
    boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
    for(size_t ch = 0; ch < channel_strings.size(); ch++){
        size_t chan = boost::lexical_cast<int>(channel_strings[ch]);
        if(chan >= usrp->get_tx_num_channels())
            throw std::runtime_error("Invalid channel(s) specified.");
        else
            channel_nums.push_back(boost::lexical_cast<int>(channel_strings[ch]));
    }
	*/


    //Lock mboard clocks
    usrp->set_clock_source(ref);

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev")) usrp->set_tx_subdev_spec(subdev);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    //set the sample rate
    if (not vm.count("rate")){
        std::cerr << "Please specify the sample rate with --rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (rate/1e6) << std::endl;
    usrp->set_tx_rate(rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...") % (usrp->get_tx_rate()/1e6) << std::endl << std::endl;

    //set the center frequency
    if (not vm.count("freq")){
        std::cerr << "Please specify the center frequency with --freq" << std::endl;
        return ~0;
    }

    for(size_t ch = 0; ch < channel_nums.size(); ch++) {
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (freq/1e6) << std::endl;
        uhd::tune_request_t tune_request(freq);
        if(vm.count("int-n")) tune_request.args = uhd::device_addr_t("mode_n=integer");
        usrp->set_tx_freq(tune_request, channel_nums[ch]);
        std::cout << boost::format("Actual TX Freq: %f MHz...") % (usrp->get_tx_freq(channel_nums[ch])/1e6) << std::endl << std::endl;

        //set the rf gain
        if (vm.count("gain")){
            std::cout << boost::format("Setting TX Gain: %f dB...") % gain << std::endl;
            usrp->set_tx_gain(gain, channel_nums[ch]);
            std::cout << boost::format("Actual TX Gain: %f dB...") % usrp->get_tx_gain(channel_nums[ch]) << std::endl << std::endl;
        }

        //set the analog frontend filter bandwidth
        if (vm.count("bw")){
            std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % bw << std::endl;
            usrp->set_tx_bandwidth(bw, channel_nums[ch]);
            std::cout << boost::format("Actual TX Bandwidth: %f MHz...") % usrp->get_tx_bandwidth(channel_nums[ch]) << std::endl << std::endl;
        }

        //set the antenna
        if (vm.count("ant")) usrp->set_tx_antenna(ant, channel_nums[ch]);
    }

    boost::this_thread::sleep(boost::posix_time::seconds(1)); //allow for some setup time

	//create a transmit streamer
	//linearly map channels (index0 = channel0, index1 = channel1, ...)
	uhd::stream_args_t stream_args("fc32");
	stream_args.channels = channel_nums;
	uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

	std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
	if (channel_nums.size() > 1) {
		// Sync times
		if (pps == "mimo") {
			UHD_ASSERT_THROW(usrp->get_num_mboards() == 2);

			//make mboard 1 a slave over the MIMO Cable
			usrp->set_time_source("mimo", 1);

			//set time on the master (mboard 0)
			usrp->set_time_now(uhd::time_spec_t(0.0), 0);

			//sleep a bit while the slave locks its time to the master
			boost::this_thread::sleep(boost::posix_time::milliseconds(100));
		} else {
			if (pps == "internal" or pps == "external" or pps == "gpsdo") {
				usrp->set_time_source(pps);
				usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
				boost::this_thread::sleep(boost::posix_time::seconds(1)); //wait for pps sync pulse
			}	
		}
	} else {
		usrp->set_time_now(0.0);
	}	
	
	//Check Ref and LO Lock detect
	std::vector<std::string> sensor_names;
	const size_t tx_sensor_chan = channel_list.empty() ? 0 : boost::lexical_cast<size_t>(channel_list[0]);
	sensor_names = usrp->get_tx_sensor_names(tx_sensor_chan);
	if (std::find(sensor_names.begin(), sensor_names.end(), "lo_locked") != sensor_names.end()) {
		uhd::sensor_value_t lo_locked = usrp->get_tx_sensor("lo_locked", tx_sensor_chan);
		std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string() << std::endl;
		UHD_ASSERT_THROW(lo_locked.to_bool());
	}
	const size_t mboard_sensor_idx = 0;
	sensor_names = usrp->get_mboard_sensor_names(mboard_sensor_idx);
	if ((ref == "mimo") and (std::find(sensor_names.begin(), sensor_names.end(), "mimo_locked") != sensor_names.end())) {
		uhd::sensor_value_t mimo_locked = usrp->get_mboard_sensor("mimo_locked", mboard_sensor_idx);
		std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string() << std::endl;
		UHD_ASSERT_THROW(mimo_locked.to_bool());
	}
	if ((ref == "external") and (std::find(sensor_names.begin(), sensor_names.end(), "ref_locked") != sensor_names.end())) {
		uhd::sensor_value_t ref_locked = usrp->get_mboard_sensor("ref_locked", mboard_sensor_idx);
		std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string() << std::endl;
		UHD_ASSERT_THROW(ref_locked.to_bool());
	}

	std::vector<int> pn_buff;
	std::vector<int> poly = {1,0,0,0,1,1,1,0,1};
	pn_buff = pn_seq_gen(poly, 255);
	/*
	std::vector<std::complex<float> > pn_buff(0);
	std::ifstream pnfile;
	std::string infilename = "PNSeq_255_MaxLenSeq.dat";
	pnfile.open(infilename.c_str(), std::ifstream::binary);
	pnfile.seekg(0, pnfile.end);
	size_t num_tx_samps = pnfile.tellg()/sizeof(std::complex<float>);
	pnfile.seekg(0, pnfile.beg);
	pn_buff.resize(num_tx_samps);
	pnfile.read((char*)&pn_buff.front(), num_tx_samps*sizeof(std::complex<float>));
	std::cout << "PN length: " << pn_buff.size() << "\n\n";
	pnfile.close();
	*/

	std::vector<std::complex<float>> pilots(fft_size - 1);
	srand(0);
	//std::cout << "Creating pilot vector...\n";
	for (int i = 0; i < fft_size - 1; i++) {
		pilots[i] = std::complex<float>((float)pn_buff[i % pn_buff.size()], 0);//(float)0.707 * std::complex<float>((float)(2*(rand() % 2) - 1), (float)(2*(rand() % 2) - 1));
	}
	std::vector<std::vector<std::complex<float>>> out_vec;

	int pre_offset = (fft_size + prefix_size)*displ[grank]*chan_per_radio;
	int post_offset = (fft_size + prefix_size)*(num_usrp*chan_per_radio - ((displ[grank] + size_of_proc_data[grank])*chan_per_radio));
	std::cout << "Rank " << grank << " Pre-offset: " << pre_offset << "\n";
	std::cout << "Rank " << grank << " Post-offset: " << post_offset << "\n";
//	if (num_usrp == 1) {
//		pre_offset = 0;
//		post_offset = 0;
//	}
	create_ofdm_sounding_frame(pilots, out_vec, fft_size, prefix_size, (int)channel_nums.size(), pre_offset, post_offset, num_reps, num_threads);
	

	std::vector<std::vector<std::complex<float> > > buff;
	buff.resize((int)channel_nums.size());
	for (int ch = 0; ch < channel_nums.size(); ch++) {
		buff[ch].resize(pn_buff.size() + out_vec[ch].size());
	}

	if (grank == 0) {
		for (int i = 0; i < pn_buff.size(); i++) {
			buff[0][i] = std::complex<float>((float)pn_buff[i], (float)pn_buff[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for (int ch = 0; ch < channel_nums.size(); ch++) {
		for (int i = 0; i < out_vec[ch].size(); i++) {
			buff[ch][pn_buff.size() + i] = out_vec[ch][i];
		}
	}


	std::vector<std::complex<float> *> buffs;
	for (int ch = 0; ch < channel_nums.size(); ch++) buffs.push_back(&buff[ch].front());
	
	std::cout << std::endl << "Tx samps: " << buff.at(0).size() << std::endl;
	
	std::signal(SIGINT, &sig_int_handler);
	std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

	uhd::tx_metadata_t md;

	md.start_of_burst = false;
	md.end_of_burst   = false;
	md.has_time_spec = true;
	std::cout << "Sending samples...\n";

	md.time_spec = usrp->get_time_now() + uhd::time_spec_t(delay);
	while(not stop_signal_called){
		//send the entire contents of the buffer
		tx_stream->send(buffs, buff.at(0).size(), md);
		md.has_time_spec = false;
	}
	md.end_of_burst = true;
	tx_stream->send("", 0, md);
	
	std::cout << "Samples sent...\n";

	//send a mini EOB packet
	//md.end_of_burst = true;
    //finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

	MPI_Finalize();
    return EXIT_SUCCESS;
}


