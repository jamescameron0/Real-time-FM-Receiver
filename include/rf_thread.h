#ifndef DY4_RF_THREAD_H
#define DY4_RF_THREAD_H

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

void rf_thread(int block_size, int block_id, std::vector<real> &block_data, std::vector<real> &i_b_half, 
        std::vector<real> &q_b_half, std::vector<real> &i_ds, std::vector<real> rf_coeff, std::vector<real> &state_i_lpf_100k,
        int rf_decim, int phase_offset_i, std::vector<real> &q_ds, std::vector<real> &state_q_lpf_100k, int phase_offset_q,
        std::vector<real> &fm_demod, float &state_I_last_sample, float &state_Q_last_sample);

#endif // DY4_RF_H