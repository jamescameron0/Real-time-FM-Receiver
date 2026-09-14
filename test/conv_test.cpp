/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

// This file shows how to write convolution unit tests, using Google C++ test framework.
// (it is based on https://github.com/google/googletest/blob/main/docs/index.md)

#include <limits.h>
#include "dy4.h"
#include "iofunc.h"
#include "filter.h"
#include "gtest/gtest.h"

namespace {

	class Convolution_Fixture: public ::testing::Test {

		public:

			const int N = 1024;	// signal size
			const int M = 101;	// tpas
			const int block_size = 128; // chunk size for block processing
			const int lower_bound = -1;
			const int upper_bound = 1;
			const real EPSILON = 1e-4;
			const int decim = 10;// decimation factor

			const int samples_to_test = 5 * block_size; // 5 blocks * blocksize is the #samples we going to use for slow method

			std::vector<real> x, h, y_reference_blk, y_ds, y_ds_sp, y_test, state_ref, state_test, y_reference_SP;
			std::vector<real> x_small;

			Convolution_Fixture() {
				x.resize(N);
				h.resize(M);
				state_ref.resize(M - 1);
				state_test.resize(M - 1);
				y_reference_blk.reserve(N+M-1);
				y_reference_SP.reserve(N+M-1);
				y_ds.reserve((N+M-1)/decim);
				y_ds_sp.reserve((N+M-1)/decim);
				y_test.reserve((N+M-1)/decim);
			}

			void SetUp() {
				generate_random_values(x, lower_bound, upper_bound);
				generate_random_values(h, lower_bound, upper_bound);

				x_small.assign(x.begin(), x.begin() + samples_to_test); //dont need to test entire data sample, 5 blocks is fien
				y_reference_blk = blockConvolve(x, h, state_ref); // full convolution using block-based implementation
				convolveFIR(y_reference_SP, x_small, h); // reference: single-pass FIR convolution (no state)
			}

			void TearDown() {
			}

			~Convolution_Fixture() {
			}
	};

	TEST_F(Convolution_Fixture, optimizedLPF_NEAR) {  //comparing regular blockConvolve + decim with optomized LPF 
		for (int i = 0; i < (int)y_reference_blk.size(); i += decim) {y_ds.push_back(y_reference_blk[i]);} //downsamplet the reference
		
		int phase_offset = 0;
		y_test = optimizedFilt(x, h, state_test, decim, phase_offset);

		for (int i = 0; i < (int)y_ds.size(); i++) {
			EXPECT_NEAR(y_ds[i], y_test[i], EPSILON) << "Original/optimizedLPF vectors differ at index " << i;
		}
		ASSERT_EQ(y_ds.size(), y_test.size()) << "Output vector sizes for decimated blockConvolve and optimizedLPF are unequal";

	}

	TEST_F(Convolution_Fixture, ployResampler_NEAR) { //comparing slow convolveFIR (single pass) with 
		y_test.clear();
		const int upS = 3;
		const int downS = 2;

		std::vector<real> x_up; // upsampled signal
		std::vector<real> y_filt; // filtered signal
		std::vector<real> y_ref; // final reference output
		bool rds = false;

		x_up.resize(x_small.size() * upS, 0.0); //increase size for upslaiming
 
		for (int i = 0; i < (int)x_small.size(); i++) { //upsample
			x_up[i * upS] = x_small[i]*upS;
		}
		convolveFIR(y_filt, x_up, h); //size = (((x * upS) + h ) - 1)/downS

		for (int i = 0; i < (int)y_filt.size(); i += downS) { //downsample
			y_ref.push_back(y_filt[i]);
		}

		std::vector<real> state_poly(h.size()/upS - 1, 0.0);
    	int phase_offset = 0;

		polyResampler(x_small, h, state_poly, upS, downS, phase_offset, y_test, rds);

		int compare_len = std::min(y_ref.size(), y_test.size()); //convolve FIR result icnludes tail, pooly doesnt, only compare the relevant data
		for (int i = 0; i < compare_len; i++) {
			EXPECT_NEAR(y_ref[i], y_test[i], EPSILON)
				<< "polyResampler differs at index " << i;
		}
		ASSERT_EQ(960, y_test.size())
       		<< "Output vector sizes mismatch in polyResampler";
	}
	TEST_F(Convolution_Fixture, optimizedLPF_BlockProcessing_OutputMatch) {

		std::vector<real> y_full;
		std::vector<real> state_full(h.size() - 1, 0.0);
		int phase_full = 0;

		y_full = optimizedFilt(x, h, state_full, decim, phase_full); // Run optimized filter on full signal

		std::vector<real> y_block; // Now process the same signal in chunks
		std::vector<real> state_block(h.size() - 1, 0.0);
		int phase_block = 0;

		for (int i = 0; i < N; i += block_size) { // Loop over signal in fixed-size blocks

			int end = std::min(i + block_size, N);
			std::vector<real> chunk(x.begin() + i, x.begin() + end);// Extract current chunk of input signal

			std::vector<real> y_chunk = optimizedFilt(  // Process chunk using same filter, preserving state + phase
				chunk, h, state_block, decim, phase_block
			);
			y_block.insert(y_block.end(), y_chunk.begin(), y_chunk.end());// Add output of block to final output
		}

		ASSERT_EQ(y_full.size(), y_block.size());

		for (int i = 0; i < (int)y_full.size(); i++) {
			EXPECT_NEAR(y_full[i], y_block[i], EPSILON)
				<< "Mismatch at index " << i;
		}
	}
	TEST_F(Convolution_Fixture, polyResampler_BlockProcessing_OutputMatch) {

		const int upS = 3; //scale factors
		const int downS = 2;
		bool rds = false;

		std::vector<real> y_full;
		std::vector<real> state_full(h.size()/upS - 1, 0.0); // Polyphase state size depends on number of phases
		int phase_full = 0;

		polyResampler(x_small, h, state_full, upS, downS, phase_full, y_full, rds);

		std::vector<real> y_block; // Process same signal in chunks 

		std::vector<real> state_block(h.size()/upS - 1, 0.0);
		int phase_block = 0;

		for (int i = 0; i < samples_to_test; i += block_size) {

			int end = std::min(i + block_size, samples_to_test);

			std::vector<real> chunk(x_small.begin() + i, x_small.begin() + end); //eah chunck 

			std::vector<real> y_chunk; //output fro thi block

			polyResampler(chunk, h, state_block, upS, downS, phase_block, y_chunk, rds); //the test ensures state and phsae are kept across blocks

			y_block.insert(y_block.end(), y_chunk.begin(), y_chunk.end()); //append output to full processing
		}

		ASSERT_EQ(y_full.size(), y_block.size());

		for (int i = 0; i < (int)y_full.size(); i++) {
			EXPECT_NEAR(y_full[i], y_block[i], EPSILON)
				<< "Mismatch at index " << i;
		}
	}

} // end of namespace