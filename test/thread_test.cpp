/*
   Comp Eng 3DY4 (Computer Systems Integration Project)
   Thread Unit Tests (GTest Version)
*/

#include <limits.h>
#include <thread>
#include <cmath>
#include "../include/dy4.h"
#include "../include/iofunc.h"
#include "../include/rf_thread.h"
#include "../include/audio_thread.h"
#include "gtest/gtest.h"

namespace {

    class Thread_Fixture: public ::testing::Test {
        public:
            const int block_size = 1024;
            const int rf_decim = 10;
            const int audio_decim = 10;
            const int audio_Fs = 48000;
            const float if_fs = 240000.0;
            const int expander = 1;
            const char format = 'm';

            std::vector<real> rf_coeff, audio_coeff, empty_vec;
            int phase_offset_i = 0, phase_offset_q = 0, phase_offset_audio = 0;

            Thread_Fixture() {
                rf_coeff.resize(101);
                audio_coeff.resize(101);
            }

            void SetUp() {
                generate_random_values(rf_coeff, -1, 1);
                generate_random_values(audio_coeff, -1, 1);
            }

            void TearDown() {}
    };

    TEST_F(Thread_Fixture, RF_Thread_Equivalence) {
        std::vector<real> block_data(block_size);
        generate_random_values(block_data, -1, 1);

        // Sequential Context
        std::vector<real> i_b_seq(block_size/2), q_b_seq(block_size/2), i_ds_seq, q_ds_seq, fm_demod_seq;
        std::vector<real> state_i_seq(100, 0.0), state_q_seq(100, 0.0);
        float state_I_last_seq = 0.0, state_Q_last_seq = 0.0;

        // Threaded Context
        std::vector<real> i_b_thr(block_size/2), q_b_thr(block_size/2), i_ds_thr, q_ds_thr, fm_demod_thr;
        std::vector<real> state_i_thr(100, 0.0), state_q_thr(100, 0.0);
        float state_I_last_thr = 0.0, state_Q_last_thr = 0.0;

        // Run Sequential
        rf_thread(block_size, 0, block_data, i_b_seq, q_b_seq, i_ds_seq, rf_coeff, state_i_seq, 
                  rf_decim, phase_offset_i, q_ds_seq, state_q_seq, phase_offset_q, 
                  fm_demod_seq, state_I_last_seq, state_Q_last_seq);

        // Run Threaded
        std::thread t_rf(rf_thread, block_size, 0, std::ref(block_data), std::ref(i_b_thr), std::ref(q_b_thr), 
                         std::ref(i_ds_thr), rf_coeff, std::ref(state_i_thr), rf_decim, phase_offset_i, 
                         std::ref(q_ds_thr), std::ref(state_q_thr), phase_offset_q, 
                         std::ref(fm_demod_thr), std::ref(state_I_last_thr), std::ref(state_Q_last_thr));
        t_rf.join();

        // Verify
        ASSERT_EQ(fm_demod_seq.size(), fm_demod_thr.size()) << "RF output sizes mismatch!";
        for (size_t i = 0; i < fm_demod_seq.size(); i++) {
            EXPECT_NEAR(fm_demod_seq[i], fm_demod_thr[i], 1e-4) << "RF Thread math mismatch at index " << i;
        }
    }

    TEST_F(Thread_Fixture, Audio_Thread_Equivalence) {
        std::vector<real> fm_demod_audio(block_size / (2 * rf_decim));
        generate_random_values(fm_demod_audio, -1, 1);

        // Sequential Context
        std::vector<real> mono_audio_seq, state_audio_seq(100, 0.0);
        std::vector<short> pcm_buf_seq(block_size / (2 * rf_decim * audio_decim));

        // Threaded Context
        std::vector<real> mono_audio_thr, state_audio_thr(100, 0.0);
        std::vector<short> pcm_buf_thr(pcm_buf_seq.size());

        // Run Sequential
        audio_thread(format, empty_vec, fm_demod_audio, empty_vec, empty_vec, 0, empty_vec, empty_vec, 
                     empty_vec, 0, empty_vec, empty_vec, empty_vec, empty_vec, empty_vec, if_fs, 
                     empty_vec, expander, mono_audio_seq, audio_coeff, state_audio_seq, audio_decim, 
                     phase_offset_audio, empty_vec, empty_vec, empty_vec, 0, 0, audio_Fs, pcm_buf_seq);

        // Run Threaded
        std::thread t_audio(audio_thread, format, std::ref(empty_vec), fm_demod_audio, empty_vec, 
                            std::ref(empty_vec), 0, std::ref(empty_vec), empty_vec, std::ref(empty_vec), 
                            0, std::ref(empty_vec), std::ref(empty_vec), std::ref(empty_vec), 
                            std::ref(empty_vec), std::ref(empty_vec), if_fs, std::ref(empty_vec), expander, 
                            std::ref(mono_audio_thr), audio_coeff, std::ref(state_audio_thr), audio_decim, 
                            phase_offset_audio, std::ref(empty_vec), empty_vec, std::ref(empty_vec), 
                            0, 0, audio_Fs, pcm_buf_thr);
        t_audio.join();

        // Verify
        ASSERT_EQ(pcm_buf_seq.size(), pcm_buf_thr.size()) << "Audio output sizes mismatch!";
        for (size_t i = 0; i < pcm_buf_seq.size(); i++) {
            EXPECT_EQ(pcm_buf_seq[i], pcm_buf_thr[i]) << "Audio Thread math mismatch at index " << i;
        }
    }

}