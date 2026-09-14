#ifndef DY4_AUDIO_THREAD_H
#define DY4_AUDIO_THREAD_H

#include "../include/dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "fm_helpers.h"

#include <thread>
#include <mutex>
#include <chrono>

void audio_thread(char format, std::vector<real> &scarFilt, std::vector<real> fm_demod, 
    std::vector<real> scar_coeff, std::vector<real> &state_scar_pilot, int scar_phase_offset,
    std::vector<real> &schFilt, std::vector<real> sch_coeff, std::vector<real> &state_sch_stereo,
    int sch_phase_offset, std::vector<real> &bufferedFm_demod, std::vector<real> &bufferedFmState,
    std::vector<real> &scarPilotI, std::vector<real> &scarPilotQ, std::vector<real> &state_scar_pll,
    float if_fs, std::vector<real> &mixedStereoFilt, int expander, std::vector<real> &mono_audio_block,
    std::vector<real> audio_coeff, std::vector<real> &state_audio_mono, int audio_decim,
    int phase_offset_audio, std::vector<real> &stereoAudioBlock, std::vector<real> stereo_audio_coeff,
    std::vector<real> &state_stereo_base, int phase_offset_stereo, int block_id, int audio_Fs, std::vector<short> pcm_buf);

#endif // DY4_AUDIO_H