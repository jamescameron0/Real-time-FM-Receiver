/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// Add headers as needed
#include "dy4.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>

struct Filters {
    std::vector<real> rf_coeff;
    std::vector<real> audio_coeff;
    std::vector<real> scar_coeff;
    std::vector<real> sch_coeff;
    std::vector<real> stereo_coeff;

    // RDS
    std::vector<real> rds_coeff;
    std::vector<real> car_rec_coeff;
    std::vector<real> demod_coeff;
    std::vector<real> rrc_coeff;
};

// INITiALIZATION FUNCTION 
void initFilters(Filters& f, float rf_Fs, float rf_Fc, int rf_taps, float if_fs, float audio_Fc, int audio_taps, float expander,
                float stereo_Fc, float scar_Fc_low, float scar_Fc_high, float sch_Fc_low, float sch_Fc_high, int rds_expander);

// IMPULSE RESPONSE FUNCTIONS
void impulseResponseLPF(real, real, unsigned short int, std::vector<real> &);
void impulseResponseBPF(real f_low, real f_high, real fs, int N_taps, std::vector<real> &h);
void impulseResponseRootRaisedCosine(const real Fs, const int N_taps, std::vector<real>& impulseResponseRRC);

// DEMOD FUNCTIONS
void fmDemod(std::vector<real> &fm_demod, const std::vector<real> &I, const std::vector<real> &Q, float &prev_I, float &prev_Q);
void fmDemodArctan(std::vector<real> &fm_demod, const std::vector<real> &I,const std::vector<real> &Q, float &prev_I, float &prev_Q);

// CONV FUNCTIONS
void convolveFIR(std::vector<real> &, const std::vector<real> &, const std::vector<real> &); //fpr testing
std::vector<real> blockConvolve(const std::vector<real> &xb, const std::vector<real> &h, std::vector<real> &state); //for testig

//efficient convolution functions 
std::vector<real> optimizedFilt(const std::vector<real> &xb, const std::vector<real> &h, std::vector<real> &state, const int decim, int &phase_offset);
void polyResampler(const std::vector<real> &xb, const std::vector<real> &h, std::vector<real> &state, const int upS, const int downS, int &phase_offset, std::vector<real> &yb, bool rds);

#endif // DY4_FILTER_H