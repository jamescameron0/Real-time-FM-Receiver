#ifndef DY4_RDS_HELPERS_H
#define DY4_RDS_HELPERS_H
#include <algorithm>
#include <cmath>

#include "../include/dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "fm_helpers.h"


void clock_data_recovery(std::vector<real> rrcWF, std::vector<real> &symbols, int sps, int &initPos, \
                        float &prevSample, float &prevSymbol, float &offset_val, float K = 0.001);
void findLocMax(std::vector<real> &rrcWindow, float &locMax, int &maxPos);


void diffDecode(std::vector<int> bits, int &state_last_bit, std::vector<int> &new_bits);

void manchesterDecode(std::vector<float> symbols, std::vector<float> &state_leftover_symbols, \
                    bool &unlocked, int &manchesterMode, float &state_last_symbol, std::vector<int> &bits);

void rds_carrier_recovery(const std::vector<real>& rds_band, float Fs, const std::vector<real>& car_rec_coeff, std::vector<real>& state_carr, 
                          std::vector<real>& state_RDS_pll, std::vector<real>& RDSPilotI, std::vector<real>& RDSPilotQ, int &rec_carr_p_off);

void rds_demodulate(const std::vector<real>& bufferedRDS_band, const std::vector<real>& RDSPilotI, const std::vector<real>& RDSPilotQ, 
                    const std::vector<real>& demod_coeff, const std::vector<real>& rrc_coeff, int expander, int decimator, std::vector<real>& state_mix, 
                    std::vector<real>& state_mixQ, std::vector<real>& state_RRC, std::vector<real>& state_RRC_Q, int& phase_offset_RRC, int& phase_offset_RRC_Q,
                    std::vector<real>& rrcFilt, std::vector<real>& rrcFiltQ, int &mixed_phase_offs, int &mixed_phase_offs_Q);

#endif