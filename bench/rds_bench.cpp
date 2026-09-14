/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

#include <benchmark/benchmark.h>
#include <thread>
#include <vector>
#include <algorithm>

#include "utils_bench.h" 
#include "../include/dy4.h"
#include "../include/rds_thread.h"
#include "../include/RDSFrameSynchro.h"
#include "../include/RDSApplicationLayer.h"

// Original block sizes from the RF bench
#define MIN_BLOCK_SIZE 96000  
#define MAX_BLOCK_SIZE 192000 

const int lower_bound = -1;
const int upper_bound = 1;

// --- Helper struct to hold all variables for the RDS thread ---
struct RDS_Context {
    int block_size;
    int rf_decim = 10; 
    
    // Core parameters
    int mode = 0;
    int block_id = 26; 
    float block_duration = 0.04f; 
    float if_FS = 240000.0f;
    int rds_expander = 1; 
    int rds_decimator = 1;

    // Offsets & Primitives
    int band_phase_offset = 0;
    int rec_carr_p_off = 0;
    int phase_offset_RRC = 0, phase_offset_RRC_Q = 0, mixed_phase_offs = 0, mixed_phase_offs_Q = 0;
    float prev_symbol = 0.0f, prev_sample = 0.0f, cdr_offset = 0.0f, state_last_symbol = 0.0f;
    int next_init_pos = 0;
    bool unlocked = true;
    int manchesterMode = 0;
    int state_last_bit = 0;

    // Vectors
    std::vector<real> fm_demod;
    
    // Filter Coefficients
    std::vector<real> rds_coeff, car_rec_coeff, demod_coeff, rrc_coeff;

    // States
    std::vector<real> state_rds, bufferedFM_state, state_carr, state_RDS_PLL;
    std::vector<real> state_mix, state_mixQ, state_RRC, state_RRC_Q;
    std::vector<real> state_leftover_symbols;
    std::vector<int> state_leftover_decoded_bits;

    // Objects
    RDSFrameSynchro rds_synchro;
    RDSApplicationLayer rds_app;

    RDS_Context(int b_size) : block_size(b_size) {
        // Size of fm_demod is block_size divided by RF decimation factor
        int fm_demod_size = block_size / (2 * rf_decim);
        fm_demod.resize(fm_demod_size, 0.0);

        // Fill input with random bench data to simulate real FM Demodulator output
        generate_random_values(fm_demod, lower_bound, upper_bound);

        // Setup filter coefficients
        rds_coeff.resize(101, 0.01);
        car_rec_coeff.resize(101, 0.01);
        demod_coeff.resize(101, 0.01);
        rrc_coeff.resize(101, 0.01);

        // Initialize States
        state_rds.resize(100, 0.0);
        bufferedFM_state.resize(50, 0.0);
        state_carr.resize(100, 0.0);
        state_RDS_PLL = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0}; // Standard PLL state size
        
        state_mix.resize(100, 0.0);
        state_mixQ.resize(100, 0.0);
        state_RRC.resize(100, 0.0);
        state_RRC_Q.resize(100, 0.0);
    }

    void reset_states() {
        std::fill(state_rds.begin(), state_rds.end(), 0.0);
        std::fill(bufferedFM_state.begin(), bufferedFM_state.end(), 0.0);
        std::fill(state_carr.begin(), state_carr.end(), 0.0);
        std::fill(state_mix.begin(), state_mix.end(), 0.0);
        std::fill(state_mixQ.begin(), state_mixQ.end(), 0.0);
        std::fill(state_RRC.begin(), state_RRC.end(), 0.0);
        std::fill(state_RRC_Q.begin(), state_RRC_Q.end(), 0.0);
        
        state_RDS_PLL = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
        state_leftover_symbols.clear();
        state_leftover_decoded_bits.clear();

        band_phase_offset = 0; rec_carr_p_off = 0;
        phase_offset_RRC = 0; phase_offset_RRC_Q = 0; mixed_phase_offs = 0;
        prev_symbol = 0.0f; prev_sample = 0.0f; cdr_offset = 0.0f; state_last_symbol = 0.0f;
        next_init_pos = 0; state_last_bit = 0; unlocked = true; manchesterMode = 0;
        
        // rds_synchro and rds_app might need reset methods depending on your implementation
        // rds_synchro.reset(); 
        // rds_app.reset();
    }
};

// --- Test 1: Measure ONLY the RDS Thread ---
static void Bench_Thread_RDS_Only(benchmark::State& state) {
    RDS_Context ctx(state.range(0));
    
    ctx.block_id = 100; 
    ctx.rds_expander = 19;
    ctx.rds_decimator = 120;
    ctx.rrc_coeff.resize(101 * ctx.rds_expander, 0.01); 

    for (auto _ : state) {
        ctx.reset_states();
        
        auto start = std::chrono::high_resolution_clock::now();
        
        rds_thread(
            ctx.mode, ctx.block_id, ctx.block_duration, ctx.fm_demod, 
            ctx.rds_coeff, ctx.state_rds, ctx.bufferedFM_state, ctx.band_phase_offset,
            ctx.if_FS, ctx.car_rec_coeff, ctx.state_carr, ctx.state_RDS_PLL, ctx.rec_carr_p_off,
            ctx.demod_coeff, ctx.rrc_coeff, ctx.state_mix, ctx.state_mixQ, 
            ctx.state_RRC, ctx.state_RRC_Q, ctx.phase_offset_RRC, ctx.phase_offset_RRC_Q, ctx.mixed_phase_offs, ctx.mixed_phase_offs_Q,
            ctx.rds_expander, ctx.rds_decimator, 
            ctx.prev_symbol, ctx.prev_sample, ctx.cdr_offset, ctx.next_init_pos, ctx.state_last_symbol, 
            ctx.state_leftover_symbols, ctx.unlocked, ctx.manchesterMode, 
            ctx.state_last_bit, 
            ctx.state_leftover_decoded_bits, ctx.rds_synchro, 
            ctx.rds_app 
        );
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        // force compiler to run the math
        benchmark::DoNotOptimize(duration_us);
        
    }
}

BENCHMARK(Bench_Thread_RDS_Only)->RangeMultiplier(2)->Range(MIN_BLOCK_SIZE, MAX_BLOCK_SIZE);