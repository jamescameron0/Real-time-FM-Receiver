#include "../include/dy4.h"
#include "../include/fm_helpers.h"

#include <fstream>

std::vector<real> delayBlock(std::vector<real> input_block, std::vector<real> &state_block){ // Delays the input signal by the length of state_block.
    std::vector<real> output_block;
    output_block.reserve(state_block.size() + (input_block.size() - state_block.size())); // output = previous state + current input (minus overlap)

    output_block.insert(output_block.end(), state_block.begin(), state_block.end());    // First part: previous samples (state)
    output_block.insert(output_block.end(), input_block.begin(), input_block.end() - state_block.size()); // Second part: current input excluding last part (which becomes new state)

    state_block.assign(input_block.end() - state_block.size(), input_block.end()); // Update state with most recent samples from input

    return output_block;
}

void fmPLL(std::vector<real> &ncoOutI, std::vector<real> &ncoOutQ, std::vector<real> &pll_states, \
            std::vector<real> pLL_in, real freq, real Fs, real ncoScale, real normBW){

    real phaseAdjust = 0.0;

    real Cp = 2.666;
    real Ci = 3.555;

    real Kp = normBW*Cp; // proportional gain
    real Ki = normBW*normBW*Ci;  // integral gain

    ncoOutI.resize(pLL_in.size() + 1);
    ncoOutQ.resize(pLL_in.size() + 1);

    real integrator = pll_states[0];
	real feedback_Phase = pll_states[1];
	real feedbackI = pll_states[2];
	real feedbackQ = pll_states[3];
	real outputPhase = pll_states[6];
    ncoOutI[0] = pll_states[4];
    ncoOutQ[0] = pll_states[5];
    
    for (size_t i = 0; i < pLL_in.size(); i++){
        real errorI, errorQ, errorD;

        errorI = pLL_in[i] * (+feedbackI); // Phase detector (I/Q mixing)
        errorQ = pLL_in[i] * (-feedbackQ);

        errorD = std::atan2(errorQ, errorI); // Computing the phase error

        integrator = integrator + Ki*errorD;

        real phase_step = 2.0 * PI * (freq/Fs) + Kp*errorD + integrator;
        feedback_Phase += phase_step;
        outputPhase += phase_step * ncoScale;

        while (feedback_Phase > 2.0 * PI) feedback_Phase -= 2.0 * PI;
        while (feedback_Phase < -2.0 * PI) feedback_Phase += 2.0 * PI;
        
        while (outputPhase > 2.0 * PI) outputPhase -= 2.0 * PI;
        while (outputPhase < -2.0 * PI) outputPhase += 2.0 * PI;

        feedbackI = std::cos(feedback_Phase);
        feedbackQ = std::sin(feedback_Phase);

        ncoOutI[i + 1] = std::cos(outputPhase + phaseAdjust);
        ncoOutQ[i + 1] = std::sin(outputPhase + phaseAdjust);  
    }
    
    pll_states = {integrator, feedback_Phase, feedbackI, feedbackQ, std::cos(outputPhase), std::sin(outputPhase), outputPhase};
}

void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<real> &block_data){
    static std::vector<unsigned char> raw_data;
    raw_data.resize(num_samples);
    std::cin.read(reinterpret_cast<char*>(raw_data.data()), num_samples);
    for (int k = 0; k < (int)num_samples; k++)
        block_data[k] = (raw_data[k] - 128) / 128.0f;
}