/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include "../include/dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "fm_helpers.h"
#include "rf_thread.h"
#include "audio_thread.h"
#include "rds_thread.h"
#include "RDSApplicationLayer.h"
#include "RDSFrameSynchro.h"

#include <thread>
#include <mutex>
#include <chrono>

int main(int argc, char** argv)
{
	float rf_Fs, audio_Fs, block_duration, expander, decimator;
	int rf_decim, audio_decim;

	float rf_Fc = 100e3;
	float  audio_Fc = 16e3;
	real stereo_Fc = 38e3;
	real scar_Fc_low = 18.5e3;
	real scar_Fc_high = 19.5e3;
	real sch_Fc_low = 22e3;
	real sch_Fc_high = 54e3;

	int mode = 0; // default
	char format = 'm'; //default
	int num_taps = 101; //constant across all modes

	//argc is the number of command line args.
	//argv is the array of strings that contains the command line args. 
    if (argc > 1) mode = atoi(argv[1]); //Check if it has at least one argument at input
	if (argc > 2) format = argv[2][0]; //Check if it has the format string, to be implemented for stereo and RDS
	if (argc > 3) num_taps = atoi(argv[3]);

	int rf_taps = num_taps; //constant across modes
	int audio_taps = num_taps;  //constant across modes
	bool rds_activate; 
	int rds_expander, rds_decimator; //expander and decimator for SPS for RDS

	//============================ MODE PARAMETER SETTINGS ======================================= // 

	switch (mode)
	{
		case  0:
			rf_Fs = 2.4e6;
			rf_decim = 10;

			audio_Fs = 48e3;
			audio_decim = 5;

			block_duration = 40e-3;

			decimator= 1;
			expander= 1;

			rds_expander = 19;
			rds_decimator = 120;
			rds_activate = (format == 'r');
		break;
		
		case  1:
			rf_Fs = 2304e3;
			rf_decim = 8;

			audio_Fs = 32e3;
			audio_decim = 9;

			block_duration = 28e-3;

			decimator= 1;
			expander= 1;

			rds_expander = 1;
			rds_decimator = 1;
			rds_activate = false;
		break;
		
		case  2:
			rf_Fs = 2.4e6;
			rf_decim = 10;

			audio_Fs = 44.1e3;
			block_duration = 60e-3;
			
			expander= 147;
			decimator= 800;
			audio_decim = decimator;

			rds_expander = 589;
			rds_decimator = 1920;
			rds_activate = (format == 'r');
		break;

		case  3:
			rf_Fs = 1152e3;
			rf_decim = 3;

			audio_Fs = 44.1e3;
			block_duration = 20e-3;
			
			expander= 147;
			decimator= 1280;
			audio_decim = decimator;

			rds_activate = false;
			rds_expander = 1;
			rds_decimator = 1;
		break;


		default: // defaults to mode 0
			mode = 0;

			rf_Fs = 2.4e6;
			rf_decim = 10;

			audio_Fs = 48e3;
			audio_decim = 5;
			
			block_duration = 40e-3;

			decimator= 1;
			expander= 1;

			rds_expander = 1;
			rds_decimator = 1;
			rds_activate = (format == 'r');
		break;
	}
	std::cerr << "Desired Block Size on mode " << mode << ": " << rf_Fs * block_duration * 2 << std::endl;

	int block_size = rf_Fs * block_duration * 2;
	float if_fs = rf_Fs/rf_decim; //intermediate frequency

	Filters filters;

	initFilters(
		filters,
		rf_Fs, rf_Fc, rf_taps,
		if_fs,
		audio_Fc, audio_taps, expander,
		stereo_Fc,
		scar_Fc_low, scar_Fc_high,
		sch_Fc_low, sch_Fc_high,
		rds_expander
	);

	std::cerr << "Operating on mode " << mode << std::endl;
	std::cerr << "Working with reals on " << sizeof(real) << " bytes" << std::endl;

	// ===================================================================================================================================================================== //
	// ================= RF, MONO, and STEREO INITIALIZATIONS  ============================================================================================================== //
	// ===================================================================================================================================================================== //
	
	std::vector<short> pcm_buf(block_size * expander / (rf_decim * audio_decim));

	std::vector<real> scarFilt, schFilt, scarPilotI, scarPilotQ, bufferedFm_demod, mixedStereoFilt;

	std::vector<real> state_i_lpf_100k(audio_taps - 1);
	std::vector<real> state_q_lpf_100k(audio_taps - 1);
	std::vector<real> state_scar_pilot(audio_taps - 1);
	std::vector<real> state_sch_stereo(audio_taps - 1);
	std::vector<real> state_stereo_base(audio_taps - 1);	
	std::vector<real> bufferedFmState((audio_taps - 1)/2);
	std::vector<real> stereoLeft;
	std::vector<real> stereoRight;

	int phase_offset_i = 0;
	int phase_offset_q = 0;
	int phase_offset_audio = 0;
	int phase_offset_stereo = 0;

	std::vector<real> state_audio_mono(audio_taps - 1);

	std::vector<real> state_scar_pll = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0};

	float state_I_last_sample = 0.0, state_Q_last_sample = 0.0;

	int sch_phase_offset = 0, scar_phase_offset = 0;


	std::vector<real> audio_data, fm_demod, mono_audio_block, stereoAudioBlock;
	std::vector<real> i_b_half(block_size/2);
	std::vector<real> q_b_half(block_size/2);

	std::vector<real> block_data(block_size);
	std::vector<real> i_ds, q_ds;

	// ===================================================================================================================================================================== //
	// ================= RDS STATE VARIABLES AND VECTORS ============================================================================== //
	// ===================================================================================================================================================================== //
	
	// CHanel Extraction
	std::vector<real> state_rds(audio_taps - 1, 0.0);
	std::vector<real> bufferedRDS_state((audio_taps - 1)/2, 0.0);
	int band_phase_offset = 0;

	// Carrier recovery
	std::vector<real> state_carr(filters.car_rec_coeff.size() - 1, 0.0);
	std::vector<real> state_RDS_PLL = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
	int rec_carr_p_off = 0;

	// Demod
	std::vector<real> state_mix(rf_taps - 1, 0.0);
	std::vector<real> state_mixQ(filters.demod_coeff.size() - 1, 0.0);

	std::vector<real> state_RRC(audio_taps* rds_expander - 1, 0.0);
	std::vector<real> state_RRC_Q(audio_taps* rds_expander - 1, 0.0);


	int phase_offset_RRC = 0;
	int phase_offset_RRC_Q = 0;
	int mixed_phase_offs = 0;
	int mixed_phase_offs_Q = 0;

	// CDR
	float prev_symbol = 0;
	float prev_sample = 0;
	float cdr_offset = 0;
	int next_init_pos = 0;
	float state_last_symbol = 0;

	// Manchester
	std::vector<real> state_leftover_symbols;
	bool unlocked = false;
	int manchesterMode = 0;
	
	// DIF DECODER
	int state_last_bit = 0;

	// Sync
	std::vector<int> state_leftover_decoded_bits;
	RDSFrameSynchro rds_synchro;

	// Application
	RDSApplicationLayer rds_app;

	// ===================================================================================================================================================================== //
	// ================= RF, AUDIO, and RDS THREAD CALLS =================================================================================================================== //
	// ===================================================================================================================================================================== //

	for (unsigned int block_id = 0; ; block_id++){	
		
		// ================= RF THREAD ====================================================== //
		std::thread rft(
			rf_thread,

			block_size, block_id,                     // params

			std::ref(block_data),                     // INPUT (raw IQ)

			std::ref(i_b_half), std::ref(q_b_half),   // IQ split 

			std::ref(i_ds), filters.rf_coeff, std::ref(state_i_lpf_100k), rf_decim, phase_offset_i, // I path

			std::ref(q_ds), std::ref(state_q_lpf_100k), phase_offset_q, // Q path

			std::ref(fm_demod),  // FM demod output

			std::ref(state_I_last_sample), // States
			std::ref(state_Q_last_sample)
		);

		rft.join();

		// ================= AUDIO THREAD ===================================================== //

		std::thread at(
			audio_thread, 
			
			format,  //mode

			std::ref(scarFilt), // SCAR - pilot filter ouput

			fm_demod, //input, fm demod

			filters.scar_coeff, std::ref(state_scar_pilot), scar_phase_offset, //scar filter coeff

    		std::ref(schFilt), filters.sch_coeff, std::ref(state_sch_stereo), sch_phase_offset, //sch (l-r band)
			
			std::ref(bufferedFm_demod), std::ref(bufferedFmState), //delayed mono path

			std::ref(scarPilotI), std::ref(scarPilotQ), //pll outputs

			std::ref(state_scar_pll), if_fs, //pll state
			
			std::ref(mixedStereoFilt), //L-R badseband
	
			expander, //for polyresampling

			std::ref(mono_audio_block), filters.audio_coeff, std::ref(state_audio_mono), audio_decim, phase_offset_audio,  // (L+R) (mono) 
			
			std::ref(stereoAudioBlock), filters.stereo_coeff, std::ref(state_stereo_base), phase_offset_stereo,  // (L-R) (stereo)
			
			block_id, audio_Fs, pcm_buf //params
		);

		std::thread rdst;

		if (rds_activate){
			rdst = std::thread(
				rds_thread,

				mode, block_id, block_duration, std::ref(fm_demod), // INPUT

				std::ref(filters.rds_coeff), std::ref(state_rds), std::ref(bufferedRDS_state), std::ref(band_phase_offset), // Channel Extraction

				if_fs, std::ref(filters.car_rec_coeff), std::ref(state_carr), std::ref(state_RDS_PLL), std::ref(rec_carr_p_off), // CARRIER RECOVERY

				std::ref(filters.demod_coeff), std::ref(filters.rrc_coeff), std::ref(state_mix), std::ref(state_mixQ), std::ref(state_RRC),// DEMODULATION

				std::ref(state_RRC_Q), std::ref(phase_offset_RRC), std::ref(phase_offset_RRC_Q), std::ref(mixed_phase_offs), std::ref(mixed_phase_offs_Q),// DEMODULATION

				rds_expander, rds_decimator, //POLYPHASE

				std::ref(prev_symbol), std::ref(prev_sample), std::ref(cdr_offset), std::ref(next_init_pos), std::ref(state_last_symbol), // CDR

				std::ref(state_leftover_symbols), std::ref(unlocked), std::ref(manchesterMode), // MANCHESTER DECODING

				std::ref(state_last_bit), // DIFF DECODER

				std::ref(state_leftover_decoded_bits), std::ref(rds_synchro), // FRAME SYNC

				std::ref(rds_app) // APP LAYER
			);
		}
		at.join();

		if (rds_activate) rdst.join(); //if we're in mode 0 or 2, rds_activate = true, we want to run them in parallel
	}
	
	/*
	cat ../data/iq_samples/samples2.raw | ./project 0 m | aplay -f S16_LE -c 1 -r 48000 -t raw
	cat ../data/iq_samples/samples0_2304000.raw | ./project 1 m | aplay -f S16_LE -c 1 -r 32000 -t raw
	cat ../data/iq_samples/samples0.raw | ./project 2 m | aplay -f S16_LE -c 1 -r 44100 -t raw
	cat ../data/iq_samples/samples0_1152000.raw | ./project 3 m | aplay -f S16_LE -c 1 -r 44100 -t raw

	cat ../data/iq_samples/stereo_l0_r9.raw | ./project 0 s | aplay -f S16_LE -c 2 -r 48000 -t raw
	cat ../data/iq_samples/stereo_l0_r9_2304000.raw | ./project 1 s | aplay -f S16_LE -c 2 -r 32000 -t raw
	cat ../data/iq_samples/stereo_l0_r9.raw | ./project 2 s | aplay -f S16_LE -c 2 -r 44100 -t raw
	cat ../data/iq_samples/stereo_l0_r9_1152000.raw | ./project 3 s | aplay -f S16_LE -c 2 -r 44100 -t raw

	cat ../data/iq_samples/samples2.raw | ./project 0 r 75 | aplay -f S16_LE -c 2 -r 48000 -t raw
	cat ../data/iq_samples/samples8.raw | ./project 2 r | aplay -f S16_LE -c 2 -r 44100 -t raw
	cat ../data/iq_samples/samples0.raw | ./project 0 r | aplay -f S16_LE -c 2 -r 48000 -t raw
	cat ../data/iq_samples/samples0.raw | ./project 2 r | aplay -f S16_LE -c 2 -r 44100 -t raw

	rtl_sdr -f 107.1M -s 2400000 - | ./project 0 m | aplay -f S16_LE -c 1 -r 48000 -t raw
	rtl_sdr -f 107.1M -s 2304000 - | ./project 1 m | aplay -f S16_LE -c 1 -r 32000 -t raw
	rtl_sdr -f 107.1M -s 2400000 - | ./project 2 m | aplay -f S16_LE -c 1 -r 44100 -t raw
	rtl_sdr -f 107.1M -s 1152000 - | ./project 3 m 83 | aplay -f S16_LE -c 1 -r 44100 -t raw

	rtl_sdr -f 107.1M -s 2400000 - | ./project 0 m | aplay -f S16_LE -c 1 -r 48000 -t raw
	rtl_sdr -f 107.1M -s 2304000 - | ./project 1 m | aplay -f S16_LE -c 1 -r 32000 -t raw
	rtl_sdr -f 107.1M -s 2400000 - | ./project 2 m | aplay -f S16_LE -c 1 -r 44100 -t raw
	rtl_sdr -f 107.1M -s 1152000 - | ./project 3 m 83 | aplay -f S16_LE -c 1 -r 44100 -t raw

	rtl_sdr -f 107.1M -s 2400000 - | ./project 0 r 75 | aplay -f S16_LE -c 2 -r 48000 -t raw
	rtl_sdr -f 107.1M -s 2304000 - | ./project 1 r 75 | aplay -f S16_LE -c 2 -r 32000 -t raw
	rtl_sdr -f 107.1M -s 2400000 - | ./project 2 r 75 | aplay -f S16_LE -c 2 -r 44100 -t raw
	rtl_sdr -f 107.1M -s 1152000 - | ./project 3 r 75 | aplay -f S16_LE -c 2 -r 44100 -t raw
	*/
	return 0;
}