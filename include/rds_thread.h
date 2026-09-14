#ifndef DY4_RDS_THREAD_H
#define DY4_RDS_THREAD_H

#include "../include/dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "fm_helpers.h"
#include "rds_helpers.h"


#include <thread>
#include <mutex>
#include <chrono>

class RDSFrameSynchro;
class RDSApplicationLayer;


void rds_thread(
    // ================= INPUT =================
    int mode,
    int block_id,
    float block_duration, 
    const std::vector<real>& fm_demod,

    // ================= CHANNEL EXTRACTION =================
    std::vector<real>& rds_coeff,
    std::vector<real>& state_rds,
    std::vector<real>& bufferedFM_state,
    int& band_phase_offset,

    // ================= CARRIER RECOVERY =================
    float if_FS,
    std::vector<real>& car_rec_coeff,
    std::vector<real>& state_carr,
    std::vector<real>& state_RDS_PLL,
    int& rec_carr_p_off,

    // ================= DEMODULATION =================
    std::vector<real>& demod_coeff,
    std::vector<real>& rrc_coeff,
    std::vector<real>& state_mix,
    std::vector<real>& state_mixQ,
    std::vector<real>& state_RRC,
    std::vector<real>& state_RRC_Q,
    int& phase_offset_RRC,
    int& phase_offset_RRC_Q,
    int& mixed_phase_offs,
    int& mixed_phase_offs_Q,

    int rds_expander, int rds_decimator, //POLY

    // ================= CDR =================
    float& prev_symbol,
    float& prev_sample,
    float& cdr_offset,
    int& next_init_pos,
    float& state_last_symbol,

    // ================= MANCHESTER =================
    std::vector<real>& state_leftover_symbols,
    bool& unlocked,
    int& manchesterMode,

    // ================= DIFF DECODER =================
    int& state_last_bit,

    // ================= FRAME SYNC =================
    std::vector<int>& state_leftover_decoded_bits,
    RDSFrameSynchro& rds_synchro,

    // ================= APPLICATION LAYER =================
    RDSApplicationLayer& rds_app
);

#endif