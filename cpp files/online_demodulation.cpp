#include <fstream>
#include <vector>
#include <algorithm>
#include <math>
#include <thread>
#include <functional>
#include <numeric>

void scale_and_shift_received_samples_for_path(std::vector<std::complex<float>> received_samples, int path, float path_gain, int path_delay, std::vector< std::vector<std::complex<float>>> &scaled_and_shifted_samples) {
	std::vector< std::complex <float> > scaled_and_shifted_samples_for_path(number_of_received_samples, std::complex<float> (0.0, 0.0));
	/* numerator of mrc -- conjugate of h times Y */
	std::transform(received_samples.begin() + path_delay, received_samples.end(), scaled_and_shifted_samples_for_path.begin(),
	               std::bind1st(std::multiplies<std::complex <float> >(), std::conj(path_gain)));
	/* denominator of mrc -- conjugate of h times h */
	std::transform(scaled_and_shifted_samples_for_path.begin(), scaled_and_shifted_samples_for_path.end(), scaled_and_shifted_samples_for_path.begin(),
	               std::bind1st(std::multiplies<float>(), 1 / (std::conj(path_gain) * path_gain)));
	scaled_and_shifted_samples[path] = scaled_and_shifted_samples_for_path;
}

int demodulate_QPSK_symbol(std::complex<float> received_constellation) {
	real_part = received_constellation.real();
	imag_part = received_constellation.imag();

	/* Computing distances of received constellation point from the reference constellation points */
	d[0] = sqrt(pow(real_part - (-1), 2) + pow(imag_part - 1, 2));
	d[1] = sqrt(pow(real_part - (-1), 2) + pow(imag_part - (-1), 2));
	d[2] = sqrt(pow(real_part - 1, 2) + pow(imag_part - 1, 2));
	d[3] = sqrt(pow(real_part - 1 , 2) + pow(imag_part - (-1), 2));

	std::vector<float> distances;
	for (int i = 0; i < 3; i++) {
		distances.push_back(d[i]);
	}
	min_distance = std::min_element(distances.begin(), distances.end());

	symbol = std::find(distances.begin(), distances.end(), min_distance) - distances.begin();

	return symbol;
}

int main() {

	int PN_SEQ_LEN = 256;
	int number_of_symbols = 1e5;
	int spread_factor = 8;
	int tail_zeros = PN_SEQ_LEN - (number_of_symbols * spread_factor) % PN_SEQ_LEN;
	int number_of_received_samples = (number_of_symbols * spread_factor) + tail_zeros;

	std::vector<float> pn_sequence; /* remember to read the PN sequence from a file and store in an array */
	std::vector<float> channel_correlations;
	std::vector<float> data_samples;

	ifstream inputFile("/received_samples.dat");

	if (inputFile) {
		double value;
		int count = 0;
		while (inputFile >> value) {
			if (count < 256) {
				channel_correlations.push_back(value);
			}
			else {
				data_samples.push_back(value);
			}
			count += 1;
		}
	}

	inputFile.close();

	std::vector< std::complex<float> > received_samples;

	for (int i = 0; i < data_samples.size(); i += 2) {
		received_samples.push_back(std::complex<float> (data_samples[i], data_samples[i + 1]));
	}

	float threshold_for_paths = 0.3;

	std::vector<float> path_gains;
	std::vector<int> path_delays;

	std::copy_if(channel_correlations.begin(), channel_correlations.end(), std::back_inserter(path_gains), std::bind(std::less<float>(), threshold_for_paths, std::placeholders::_1));

	for (unsigned i = 0; i < path_gains.size(); i++) {
		path_delays.push_back(std::find(channel_correlations.begin(), channel_correlations.end(), path_gains[i]) - channel_correlations.begin());
	}

	/* Assume that data stored in file is of the format real, imaginary, real, imaginary, and so on */

	/* For each path, multiply the received samples by the complex channel gain and shift by delay */
	std::vector<std::vector<std::complex<float>>> scaled_and_shifted_samples(path_gains.size(), std::vector<std::complex <float> >(number_of_received_samples, std::complex<float> (0.0, 0.0)));
	std::thread threads[path_gains.size()];
	for (int path = 0; path < path_gains.size(); path++) {
		threads[path] = thread(scale_and_shift_received_samples_for_path, received_samples, path, path_gains[path], path_delays[path], scaled_and_shifted_samples);
	}
	for (auto& th : threads) {
		th.join();
	}

	/* Computing the sum over all paths to yield a vector of size PN_SEQ_LEN */
	std::vector<std::complex<float>> summed_samples(number_of_received_samples, 0.0);
	int i = 0;
	while (i < scaled_and_shifted_samples.size()) {
		std::transform(summed_samples[i].begin(), summed_samples[i].end(), scaled_and_shifted_samples[i].begin(),
		               summed_samples.begin(), std::plus<float>());

	}

	/* Construct vector of PN sequence (replicated) of length = number_of_received_samples */

	auto number_of_repetitions = number_of_received_samples / PN_SEQ_LEN;
	std::vector<std::complex<float>> pn_sequence_replicated(number_of_received_samples);
	for (int i = 0; i < number_of_repetitions; i++) {
		std::copy(pn_sequence.begin(), pn_sequence.end(), pn_sequence_replicated.begin() + (i * PN_SEQ_LEN));
	}


	/* Despreading by computing multiplying the elements of PN sequence and the scaled + shifted received samples */

	std::vector<std::complex<float>> despreaded_samples(number_of_received_samples, std::complex<float> (0.0, 0.0));

	std::transform(scaled_and_shifted_samples.begin(), scaled_and_shifted_samples.end(), pn_sequence_replicated.begin(), despreaded_samples.begin(), std::multiplies<std::complex<float>>() );

	/* Delete the tail zeros (added because of spread factor) */

	despreaded_samples.erase(despreaded_samples.end(), despreaded_samples.end() - tail_zeros);

	/* Average 'spread_factor' number of samples to obtain a number_of_symbols x 1 vector */

	std::vector<std::complex<float>> received_constellation(number_of_symbols, std::complex<float>(0.0, 0.0));

	for (int i = 0; i < number_of_symbols; i++) {
		received_constellation[i] = std::accummulate(despreaded_samples.begin() + (i * spread_factor), despreaded_samples.begin() + (i * spread_factor) + spread_factor, std::complex<float>(0.0, 0.0));
	}

	std::transform(received_constellation.begin(), received_constellation.end(), received_constellation.begin(), std::bind1st(std::multiplies<float>(), 1 / spread_factor));

	std::vector <int> demodulated_QPSK_symbols (number_of_symbols);
	std::transform(received_constellation.begin(), received_constellation.end(), demodulated_QPSK_symbols.begin(), [](std::complex<float> value) {return demodulate_QPSK_symbol(value)});
	return 0;
}
