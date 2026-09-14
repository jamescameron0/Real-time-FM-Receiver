#include "../include/dy4.h"
#include "../include/rds_helpers.h"
#include "filter.h"

#include <iostream>
#include <cmath>

void rds_carrier_recovery(const std::vector<real>& rds_band, float Fs, const std::vector<real>& car_rec_coeff, std::vector<real>& state_carr,
                          std::vector<real>& state_RDS_pll, std::vector<real>& RDSPilotI, std::vector<real>& RDSPilotQ, int &rec_carr_p_off) {

    // 1. Squaring Non-Linearity to double the phase and remove modulation
    std::vector<real> rds_band_squared(rds_band.size());
    for (size_t i = 0; i < rds_band.size(); i++) {
        rds_band_squared[i] = rds_band[i] * rds_band[i];
    }

    // 2. Bandpass Filter (114 kHz) using blockConvolve
    std::vector<real> rds_recovered_carr = optimizedFilt(rds_band_squared, car_rec_coeff, state_carr, 1, rec_carr_p_off);

    // 3. PLL with NCO to generate the 57 kHz pilot tones
    fmPLL(RDSPilotI, RDSPilotQ, state_RDS_pll, rds_recovered_carr, 114e3, Fs, 0.5, 0.001); 
}

void rds_demodulate(const std::vector<real>& bufferedRDS_band, const std::vector<real>& RDSPilotI, const std::vector<real>& RDSPilotQ, 
                    const std::vector<real>& demod_coeff, const std::vector<real>& rrc_coeff, int expander, int decimator, 
                    std::vector<real>& state_mix, std::vector<real>& state_mixQ, std::vector<real>& state_RRC, std::vector<real>& state_RRC_Q, 
                    int& phase_offset_RRC, int& phase_offset_RRC_Q,std::vector<real>& rrcFilt, std::vector<real>& rrcFiltQ, int &mixed_phase_offs, int &mixed_phase_offs_Q) {

    size_t block_len = std::min(bufferedRDS_band.size(), RDSPilotI.size());
    bool rds_resample = true;
    // 1. Mix to baseband
    std::vector<real> mixedRDSsignalI(block_len);
    std::vector<real> mixedRDSsignalQ(block_len);
    for (size_t i = 0; i < block_len; i++) {
        mixedRDSsignalI[i] = RDSPilotI[i] * bufferedRDS_band[i];
        mixedRDSsignalQ[i] = RDSPilotQ[i] * bufferedRDS_band[i];
    }

    // 2. Lowpass Filtering using optFilt
    std::vector<real> mixedFilt = optimizedFilt(mixedRDSsignalI, demod_coeff, state_mix, 1, mixed_phase_offs);
    std::vector<real> mixedFiltQ = optimizedFilt(mixedRDSsignalQ, demod_coeff, state_mixQ, 1, mixed_phase_offs_Q);

    // 3. Rational Resampling (Upsample -> RRC Filter -> Decimate)
    polyResampler(mixedFilt, rrc_coeff, state_RRC, expander, decimator, phase_offset_RRC, rrcFilt, rds_resample);   //240k -> 16 x 2375 or 31 x 2375 Sample rate
    polyResampler(mixedFiltQ, rrc_coeff, state_RRC_Q, expander, decimator, phase_offset_RRC_Q, rrcFiltQ, rds_resample);
}

void clock_data_recovery(std::vector<real> rrcWF, std::vector<real> &symbols, int sps, int &initPos, float &prevSample, float &prevSymbol, \
                        float &offset_val, float K){
    float offset = offset_val;
    int i = initPos;
    float curr_sample, curr_symbol, timing_err;

    symbols.clear();

    while (i < (int)rrcWF.size())
    {
        curr_sample = rrcWF[i];
        curr_symbol = (curr_sample > 0) ? 1.0 : -1.0;

        if (prevSymbol != curr_symbol) timing_err = prevSymbol * curr_sample - curr_symbol * prevSample;
        else timing_err = 0.0;

        offset += timing_err * K;
        offset = std::clamp(offset, (float)(-sps/2), (float)(sps/2));

        symbols.push_back(curr_sample);

        prevSample = curr_sample;
        prevSymbol = curr_symbol;

        i += sps + static_cast<int>(std::round(offset));
    }
    initPos = i - (int)rrcWF.size();
    offset_val = offset;
}


void findLocMax(std::vector<real> &rrcWindow, float &locMax, int &maxPos){
    locMax = 0;
    float pt;

    for (size_t i = 0; i < rrcWindow.size(); i++)
    {
        pt = std::abs(rrcWindow[i]);
        if (pt > locMax)
        {
            locMax = pt;
            maxPos = i;
        }
    }
}

void diffDecode(std::vector<int> bits, int &state_last_bit, std::vector<int> &new_bits){ //current bits, last bit of prev block, decoded bits
    new_bits.clear();

    if (bits.empty()) return; //if no input nothging to decode

    new_bits.push_back(bits[0] ^ state_last_bit); //just put bits ^ bits-1=last_bit in new_bits first

    for (size_t i = 1; i < bits.size(); i++)
    {
        new_bits.push_back(bits[i] ^ bits[i - 1]);
    }
    state_last_bit = bits.back(); //bits[-1]; this would be in python .back returns last element
 }







// Must adjust this to support relative HL detection
void manchesterDecode(std::vector<float> symbols, std::vector<float> &state_leftover_symbols, \
                    bool &unlocked, int &manchesterMode, float &state_last_symbol, std::vector<int> &bits){
    std::vector<float> allSymbols;
    std::vector<int> bits0, bits1;
    bits.clear();
    
    const float threshold = 0.25; // Configurable threshold value to determine whether it is a transition

    if (!std::isnan(state_last_symbol)) //nan is special value called Not-a-Number (NaN). isnan returns if tgis value is NAN
    {
        allSymbols.push_back(state_last_symbol); //cant assign to index 0 since it is empty
        allSymbols.insert(allSymbols.end(), symbols.begin(), symbols.end());
    }
    else {
        if ((int)state_leftover_symbols.size() > 0)
        {
            allSymbols = state_leftover_symbols;
            allSymbols.insert(allSymbols.end(), symbols.begin(), symbols.end());
        }
        else allSymbols = symbols;
    }
    
    state_leftover_symbols.clear();
    int pos = 0, desync = 0;
    
    if (unlocked)
    {
        int desync0 = 0, desync1 = 0, k = 0;
        
        // 0 offset
        while (k < (int)allSymbols.size() - 1 && k < 24)
        {
            float s0 = allSymbols[k];
            float s1 = allSymbols[k + 1];
            float diff = s0 - s1;
            float local_amp = std::abs(s0) + std::abs(s1);
            float normalized_diff = (local_amp > 1e-6f) ? (diff / local_amp) : 0.0f;
            
            if (normalized_diff > threshold) // HL 
            {
                bits0.push_back(1);
            }
            else if (normalized_diff < -threshold) // LH
            {
                bits0.push_back(0);
            }
            else
            {
                desync0++;
            }
            k += 2;        
        }

        // 1 offset
        k = 1;
        while (k < (int)allSymbols.size() - 1 && k < 23)
        {
            float s0 = allSymbols[k];
            float s1 = allSymbols[k + 1];
            float diff = s0 - s1;
            float local_amp = std::abs(s0) + std::abs(s1);
            float normalized_diff = (local_amp > 1e-6f) ? (diff / local_amp) : 0.0f;
            
            if (normalized_diff > threshold)
            {
                bits1.push_back(1);
            }
            else if (normalized_diff < -threshold)
            {
                bits1.push_back(0);
            }
            else
            {
                desync1++;
            }
            k += 2;        
        }
        // Check which one had less desyncs (HH or LL_) and set the 'manchesterMode' accordingly
        manchesterMode = (desync0 <= desync1) ? 0 : 1;
        if (manchesterMode == 0) bits.insert(bits.end(), bits0.begin(), bits0.end());
        else if (manchesterMode == 1) bits.insert(bits.end(), bits1.begin(), bits1.end());

        pos = 24 - manchesterMode;
        unlocked = false;
    }
    
    while (pos < (int)allSymbols.size() - 1)
    {
        float s0 = allSymbols[pos];
        float s1 = allSymbols[pos + 1];
        float diff = s0 - s1;
        float local_amp = std::abs(s0) + std::abs(s1);
        float normalized_diff = (local_amp > 1e-6f) ? (diff / local_amp) : 0.0f;

        if (normalized_diff > threshold) // HL transition
        {
            bits.push_back(1);
            pos += 2;
            desync = 0;
        }
        else if (normalized_diff < -threshold) // LH transition
        {
            bits.push_back(0);
            pos += 2;
            desync = 0;
        }
        else
        {
            desync++;
            pos += 2;
            if (desync > 1)
            {
                unlocked = true;
                state_last_symbol = NAN;
                int safe_idx = std::max(0, pos - 5);
                state_leftover_symbols.insert(state_leftover_symbols.begin(), allSymbols.begin() + safe_idx, allSymbols.end());
                break;
            }
            
        } 
    }
    if (pos == (int)allSymbols.size() - 1) state_last_symbol = allSymbols.back();
    else state_last_symbol = NAN;
}