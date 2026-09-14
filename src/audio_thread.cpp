#include "../include/dy4.h"
#include "../include/audio_thread.h"

void audio_thread(char format, std::vector<real> &scarFilt, std::vector<real> fm_demod, 
    std::vector<real> scar_coeff, std::vector<real> &state_scar_pilot, int scar_phase_offset,
    std::vector<real> &schFilt, std::vector<real> sch_coeff, std::vector<real> &state_sch_stereo,
    int sch_phase_offset, std::vector<real> &bufferedFm_demod, std::vector<real> &bufferedFmState,
    std::vector<real> &scarPilotI, std::vector<real> &scarPilotQ, std::vector<real> &state_scar_pll,
    float if_fs, std::vector<real> &mixedStereoFilt, int expander, std::vector<real> &mono_audio_block,
    std::vector<real> audio_coeff, std::vector<real> &state_audio_mono, int audio_decim,
    int phase_offset_audio, std::vector<real> &stereoAudioBlock, std::vector<real> stereo_audio_coeff,
    std::vector<real> &state_stereo_base, int phase_offset_stereo, int block_id, int audio_Fs, std::vector<short> pcm_buf){
    
    // STEREO PRCOESSING
    if (format == 's' || format == 'r'){ 
        scarFilt = optimizedFilt(fm_demod, scar_coeff, state_scar_pilot, 1, scar_phase_offset); // Extract the 19 kHz pilot tone (SCAR band)

        schFilt = optimizedFilt(fm_demod, sch_coeff, state_sch_stereo, 1, sch_phase_offset); // Extract the stereo diffrence signal (L - R) centered at 38 kHz (SCH band)

        bufferedFm_demod = delayBlock(fm_demod, bufferedFmState);  // Delay the FM demodulated signal to align timing with stereo path

        fmPLL(scarPilotI, scarPilotQ, state_scar_pll, scarFilt, 19e3, if_fs, 2.0, 0.01); // Run PLL on the 19 kHz polot to recover a clean carrier (arojnd 38kHz)

        size_t mix_len = std::min(scarPilotI.size() - 1, schFilt.size()); // Mix recovered carrier with (L - R) signal
        mixedStereoFilt.clear();
        for (size_t i = 0; i < mix_len; i++) mixedStereoFilt.push_back(2 * scarPilotI[i] * schFilt[i]);
    }

    bool rds_resample = false;

    // ADUIO FILTERING //

    if (expander == 1) { //decimation depending on expander factor

        if (format == 'm')
        {
            mono_audio_block = optimizedFilt(fm_demod, audio_coeff, state_audio_mono, audio_decim, phase_offset_audio); //mono extract L+R
        }
        else if (format == 's' || format == 'r'){

            // These two are independent once bufferedFm_demod and mixedStereoFilt are ready, so run in parallel 
            std::thread mono_t([&](){
                mono_audio_block = optimizedFilt(bufferedFm_demod, audio_coeff, state_audio_mono, audio_decim, phase_offset_audio); // FILTER mono
            });
            std::thread stereo_t([&](){
                stereoAudioBlock = optimizedFilt(mixedStereoFilt, stereo_audio_coeff, state_stereo_base, audio_decim, phase_offset_stereo); // Filter baseband (L - R)
            });

            mono_t.join();
            stereo_t.join();
        }
    }
    else { // If mode 2 or 3, apply resampling to decimate to 44.1 kHz
    
        if (format == 'm')
            polyResampler(fm_demod, audio_coeff, state_audio_mono, expander, audio_decim, phase_offset_audio, mono_audio_block, rds_resample);
        
        else if (format == 's' || format == 'r'){
            // These two are independent once bufferedFm_demod and mixedStereoFilt are ready, so run in parallel
            std::thread mono_t([&](){
                polyResampler(bufferedFm_demod, audio_coeff, state_audio_mono, 
                            expander, audio_decim, phase_offset_audio, mono_audio_block, rds_resample); 
            });                                                                                 // same as before but with poly
            std::thread stereo_t([&](){
                polyResampler(mixedStereoFilt, stereo_audio_coeff, state_stereo_base, 
                            expander, audio_decim, phase_offset_stereo, stereoAudioBlock, rds_resample);
            });

            mono_t.join();
            stereo_t.join();
        }
    }

    // ================= OUTPUT =================

    if (format == 'm'){ 
        for (const auto& sample : mono_audio_block) {  // Mono output: write (L + R)
            short pcm_sample = static_cast<short>(std::max(-1.0f, std::min(1.0f, (float)sample)) * 16384.0f);
            //comment out this line when running benchmark   
            fwrite(&pcm_sample, sizeof(short), 1, stdout);
        }
    }
    else if (format == 's' || format == 'r'){
        
        size_t len = std::min(mono_audio_block.size(), stereoAudioBlock.size()); 
                                           // Combine (L+R) and (L-R) to recover Left and Right channels:
        for (size_t i = 0; i < len; i++) { // (L+R) + (L-R) = 2L x 0.5 = L // (L+R) - (L-R) = 2R x 0.5 = R
            pcm_buf[2*i]   = static_cast<short>(std::max(-1.0f, std::min(1.0f, (float)(0.5 * (mono_audio_block[i] + stereoAudioBlock[i]))))  * 16384.0f); // LEFT
            pcm_buf[2*i+1] = static_cast<short>(std::max(-1.0f, std::min(1.0f, (float)(0.5 * (mono_audio_block[i] - stereoAudioBlock[i])))) * 16384.0f); // RIGHT
        }
            //comment out these 2 lines when running benchmark
           fwrite(pcm_buf.data(), sizeof(short), pcm_buf.size(), stdout);
           fflush(stdout);
    }
}