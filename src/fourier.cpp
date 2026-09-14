/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// Source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

// Just DFT function (no FFT)
void DFT(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * PI * (k * m) / x.size());
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// Function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<real>> &Xf, std::vector<real> &Xmag)
{
	Xmag.clear();
	Xmag.resize(Xf.size(), real(0));
	for (int i = 0; i < (int)Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i]) / Xf.size();
	}
}
// Add your own code to estimate the PSD

void estimatePSD(const std::vector<real> &samples, real Fs, std::vector<real> &freq, std::vector<real> &psd_est) {
	// rename the NFFT argument (notation consistent with matplotlib.psd)
	// to freq_bins (i.e., frequency bins for which we compute the spectrum)
	int freq_bins = NFFT;
	// frequency increment (or resolution of the frequency bins)
	real df = Fs / freq_bins;

	// create the frequency vector to be used on the X axis
	// for plotting the PSD on the Y axis (only positive freq)
	freq.clear();
	for (int i = 0; i < Fs / 2; ++i) {
		freq.push_back(i * df);
	}

	// design the Hann window used to smoothen the discrete data in order
	// to reduce the spectral leakage after the Fourier transform
	std::vector<real> hann(freq_bins);
	for (int i = 0; i < freq_bins; ++i) {
		hann[i] = 0.5 * (1 - std::cos(2 * PI * i / (freq_bins - 1)));
	}

	// samples should be a multiple of frequency bins, so
	// the number of segments used for estimation is an integer
	// note: for this to work you must provide an argument for the
	// number of frequency bins not greater than the number of samples!		
	int no_segments = static_cast<int>(std::floor(samples.size() / static_cast<real>(freq_bins)));
	if (no_segments < 1) {
		psd_est.clear();
		return;
	}

	// create an empty list where the PSD for each segment is computed
	std::vector<real> psd_list;

	// iterate through all the segments
	for (int k = 0; k < no_segments; ++k) {

		// apply the hann window (using pointwise multiplication)
		// before computing the Fourier transform on a segment
		std::vector<real> windowed_samples(freq_bins);
		for (int i = 0; i < freq_bins; ++i) {
			windowed_samples[i] = samples[k * freq_bins + i] * hann[i];
		}

		// since input is real, we keep only the positive half of the spectrum
		// however, we will also add the signal energy of negative frequencies
		// to have a better and more accurate PSD estimate when plotting
		std::vector<std::complex<real>> Xf;
		DFT(windowed_samples, Xf);

		// keep only positive freq bins
		Xf.resize(freq_bins / 2);

		std::vector<double> psd_seg(Xf.size());
		for (size_t i = 0; i < Xf.size(); ++i) {
			// compute signal power
			// add the energy from the negative freq bins
			psd_seg[i] = 2.0 * 1.0 / (Fs * freq_bins / 2.0) * std::norm(Xf[i]);
		}

		// append to the list where PSD for each segment is stored
		// in sequential order (first segment, followed by the second one, ...)
		psd_list.insert(psd_list.end(), psd_seg.begin(), psd_seg.end());
	}

	// iterate through all the frequency bins (positive freq only)
	// from all segments and average them (one bin at a time ...)
	std::vector<real> psd_avg(freq_bins / 2, 0.0);
	for (int k = 0; k < freq_bins / 2; ++k) {
		// iterate through all the segments
		for (int l = 0; l < no_segments; ++l) {
			psd_avg[k] += psd_list[k + l * static_cast<int>(freq_bins/2)];
		}
		// compute the estimate for each bin
		psd_avg[k] /= no_segments;
	}

	// translate to the decibel (dB) scale
	psd_est.resize(freq_bins / 2);
	for (int k = 0; k < freq_bins / 2; ++k) {
		psd_est[k] = 10.0f * std::log10(psd_avg[k] + 1e-20f);
	}
}

void DFT_reference(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / x.size());
			Xf[m] +=  + x[k] * std::exp(expval);
		}
	}
}

void DFT_init_bins(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / N);
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// ALGO 2 //
void DFT_code_motion(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size(); 
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	std::complex<real> base_exponent(0, -2 * M_PI / N);
	for (int m = 0; m < N; m++) {
		std::complex<real> intermediate_exponent = base_exponent * (real)m;
		for (int k = 0; k < N; k++) {
			std::complex<real> twiddle_exponent = intermediate_exponent * (real)k;
			Xf[m] += x[k] * std::exp(twiddle_exponent);
		}
	}
}

// ALGO 3 //

void generate_DFT_twiddles(const int& N, std::vector<std::complex<real>> &Twiddle1D) {

	Twiddle1D.resize(N);
	for (int k = 0; k < N; k++) {
		std::complex<real> expval(0, -2 * M_PI * k / N);
		Twiddle1D[k] = std::exp(expval);
	}
}

void DFT_precomp_twiddle(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();
	std::vector<std::complex<real>> Twiddle1D;

	generate_DFT_twiddles(N, Twiddle1D);
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			Xf[m] += x[k] * Twiddle1D[(k*m) % N];
		}
	}
}

// ALGO 4 //

void generate_DFT_matrix(const int& N, std::vector<std::vector<std::complex<real>>> &Twiddle2D) {

	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));
    std::vector<std::complex<real>> Twiddle1D;
	generate_DFT_twiddles(N, Twiddle1D);

	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			Twiddle2D[m][k] = Twiddle1D[(k * m) % N];
		}
	}
}

void DFT_precomp_matrix(const std::vector<real> &x, std::vector<std::complex<real>> &Xf, const std::vector<std::vector<std::complex<real>>> &Twiddle2D){
	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	for (int m = 0; m < N; m++) {
		for(int k =0; k<N; k++){
			Xf[m] += x[k] * Twiddle2D[m][k];
		}
	}
}

// ALGO 5 ///

void DFT_precomp_matrix_reordered(const std::vector<real> &x, std::vector<std::complex<real>> &Xf, const std::vector<std::vector<std::complex<real>>> &Twiddle2D){
	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	for (int k = 0; k < N; k++) {
		for(int m =0; m<N; m++){
			Xf[m] += x[k] * Twiddle2D[m][k];
		}
	}
}

// ALGO 6 ///

void DFT_precomp_loop_unrolling(const std::vector<real> &x, std::vector<std::complex<real>> &Xf, const std::vector<std::vector<std::complex<real>>> &Twiddle2D){

	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	int unrolled_limit = N - (N%4);
	for (int m = 0; m < N; m++) {
		for(int k =0; k<unrolled_limit; k+=4){
			Xf[m] += x[k] * Twiddle2D[m][k] +  x[k+1] * Twiddle2D[m][k+1] +  x[k+2] * Twiddle2D[m][k+2] +  x[k+3] * Twiddle2D[m][k+3];
		}
		for(int k = unrolled_limit; k<N; k++){
			Xf[m] += x[k] * Twiddle2D[m][k];
		}
	}

}