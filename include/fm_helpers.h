
#ifndef DY4_FM_HELPERS_H
#define DY4_FM_HELPERS_H

// Add headers as needed
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>


std::vector<real> delayBlock(std::vector<real> input_block, std::vector<real> &state_block);

void fmPLL(std::vector<real> &ncoOutI, std::vector<real> &ncoOutQ, std::vector<real> &pll_states, \
            std::vector<real> pLL_in, real freq, real Fs, real ncoScale, real normBW);
            
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<real> &block_data);
#endif