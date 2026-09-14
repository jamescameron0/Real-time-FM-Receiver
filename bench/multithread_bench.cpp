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
#include <chrono>

#include "utils_bench.h" 
#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "rf_thread.h"
#include "audio_thread.h"

// We will test block sizes representing Mode 0
#define MIN_BLOCK_SIZE 96000  // 2.4M * 20ms * 2
#define MAX_BLOCK_SIZE 192000 // 2.4M * 40ms * 2

const int lower_bound = -1;
const int upper_bound = 1;

// --- Helper struct to hold all the variables from your main.cpp ---
struct SDR_Context {
    int block_size;
    int rf_decim = 10;
    int audio_decim = 5;
    float expander = 1.0;
    float if_fs = 240000.0;
    float audio_Fs = 48000.0;
    int block_id = 0;
    char format = 's'; // Testing stereo 's'

    int phase_offset_i = 0, phase_offset_q = 0, phase_offset_audio = 0, phase_offset_stereo = 0;
    int sch_phase_offset = 0, scar_phase_offset = 0;
    float state_I_last_sample = 0.0, state_Q_last_sample = 0.0;

    std::vector<real> block_data, i_b_half, q_b_half, i_ds, q_ds;
    std::vector<real> fm_demod, fm_demod_pipeline_in;

    std::vector<real> rf_coeff, audio_coeff, scar_coeff, sch_coeff, stereo_audio_coeff;
    std::vector<real> state_i_lpf_100k, state_q_lpf_100k, state_scar_pilot, state_sch_stereo;
    std::vector<real> state_stereo_base, state_audio_mono, bufferedFmState, state_scar_pll;

    std::vector<real> scarFilt, schFilt, bufferedFm_demod, scarPilotI, scarPilotQ;
    std::vector<real> mixedStereoFilt, mono_audio_block, stereoAudioBlock;
    std::vector<short> pcm_buf;

    SDR_Context(int b_size) : block_size(b_size) {
        // Allocate raw data buffers
        block_data.resize(block_size, 0.0);
        i_b_half.resize(block_size / 2, 0.0);
        q_b_half.resize(block_size / 2, 0.0);
        fm_demod.resize(block_size / (2 * rf_decim), 0.0);
        fm_demod_pipeline_in.resize(block_size / (2 * rf_decim), 0.0);

        // Fill input with random bench data
        generate_random_values(block_data, lower_bound, upper_bound);
        generate_random_values(fm_demod_pipeline_in, lower_bound, upper_bound);

        // Setup filter coefficients 
        rf_coeff.resize(101, 0.01);
        audio_coeff.resize(101, 0.01);
        scar_coeff.resize(101, 0.01);
        sch_coeff.resize(101, 0.01);
        stereo_audio_coeff.resize(101, 0.01);

        // Setup States
        state_i_lpf_100k.resize(100, 0.0);
        state_q_lpf_100k.resize(100, 0.0);
        state_scar_pilot.resize(100, 0.0);
        state_sch_stereo.resize(100, 0.0);
        state_stereo_base.resize(100, 0.0);
        state_audio_mono.resize(100, 0.0);
        bufferedFmState.resize(50, 0.0);
        state_scar_pll = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
        
        pcm_buf.resize(block_size * expander / (rf_decim * audio_decim), 0);
    }

    void reset_states() {
        std::fill(state_i_lpf_100k.begin(), state_i_lpf_100k.end(), 0.0);
        std::fill(state_q_lpf_100k.begin(), state_q_lpf_100k.end(), 0.0);
        std::fill(state_scar_pilot.begin(), state_scar_pilot.end(), 0.0);
        std::fill(state_sch_stereo.begin(), state_sch_stereo.end(), 0.0);
        std::fill(state_stereo_base.begin(), state_stereo_base.end(), 0.0);
        std::fill(state_audio_mono.begin(), state_audio_mono.end(), 0.0);
        std::fill(bufferedFmState.begin(), bufferedFmState.end(), 0.0);
        state_scar_pll = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
        
        phase_offset_i = 0; phase_offset_q = 0; 
        phase_offset_audio = 0; phase_offset_stereo = 0;
        sch_phase_offset = 0; scar_phase_offset = 0;
        state_I_last_sample = 0.0; state_Q_last_sample = 0.0;
    }
};

// --- Test 1: Measure ONLY the RF Thread ---
static void Bench_Thread_RF_Only(benchmark::State& state) {
    SDR_Context ctx(state.range(0));
    ctx.block_id = 100; // Bypass startup locks
    
    for (auto _ : state) {
        ctx.reset_states();
        
        auto start = std::chrono::high_resolution_clock::now();

        rf_thread(ctx.block_size, ctx.block_id, ctx.block_data, ctx.i_b_half, 
            ctx.q_b_half, ctx.i_ds, ctx.rf_coeff, ctx.state_i_lpf_100k,
            ctx.rf_decim, ctx.phase_offset_i, ctx.q_ds, ctx.state_q_lpf_100k, ctx.phase_offset_q,
            ctx.fm_demod, ctx.state_I_last_sample, ctx.state_Q_last_sample);
            
        auto end = std::chrono::high_resolution_clock::now();
        auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        benchmark::DoNotOptimize(duration_us);

        benchmark::ClobberMemory(); 
    }
}
BENCHMARK(Bench_Thread_RF_Only)->RangeMultiplier(2)->Range(MIN_BLOCK_SIZE, MAX_BLOCK_SIZE);

// --- Test 2: Measure ONLY the Audio Thread ---
static void Bench_Thread_Audio_Only(benchmark::State& state) {
    SDR_Context ctx(state.range(0));
    ctx.block_id = 100;
    
    for (auto _ : state) {
        ctx.reset_states();
        
        auto start = std::chrono::high_resolution_clock::now();

        audio_thread(ctx.format, ctx.scarFilt, ctx.fm_demod_pipeline_in, ctx.scar_coeff, ctx.state_scar_pilot, ctx.scar_phase_offset,
            ctx.schFilt, ctx.sch_coeff, ctx.state_sch_stereo, ctx.sch_phase_offset, ctx.bufferedFm_demod, ctx.bufferedFmState, ctx.scarPilotI, ctx.scarPilotQ,
            ctx.state_scar_pll, ctx.if_fs, ctx.mixedStereoFilt, ctx.expander, ctx.mono_audio_block, ctx.audio_coeff, ctx.state_audio_mono, ctx.audio_decim,
            ctx.phase_offset_audio, ctx.stereoAudioBlock, ctx.stereo_audio_coeff, ctx.state_stereo_base, ctx.phase_offset_stereo, ctx.block_id, ctx.audio_Fs, ctx.pcm_buf);
            
        auto end = std::chrono::high_resolution_clock::now();
        auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        benchmark::DoNotOptimize(duration_us);

        benchmark::ClobberMemory(); 
    }
}
BENCHMARK(Bench_Thread_Audio_Only)->RangeMultiplier(2)->Range(MIN_BLOCK_SIZE, MAX_BLOCK_SIZE);

// --- Test 3: Sequential ---
static void Bench_Sequential(benchmark::State& state) {
    SDR_Context ctx(state.range(0));
    ctx.block_id = 100;
    
    for (auto _ : state) {
        ctx.reset_states();
        
        auto start = std::chrono::high_resolution_clock::now();

        rf_thread(ctx.block_size, ctx.block_id, ctx.block_data, ctx.i_b_half, 
            ctx.q_b_half, ctx.i_ds, ctx.rf_coeff, ctx.state_i_lpf_100k,
            ctx.rf_decim, ctx.phase_offset_i, ctx.q_ds, ctx.state_q_lpf_100k, ctx.phase_offset_q,
            ctx.fm_demod, ctx.state_I_last_sample, ctx.state_Q_last_sample);

        audio_thread(ctx.format, ctx.scarFilt, ctx.fm_demod, ctx.scar_coeff, ctx.state_scar_pilot, ctx.scar_phase_offset,
            ctx.schFilt, ctx.sch_coeff, ctx.state_sch_stereo, ctx.sch_phase_offset, ctx.bufferedFm_demod, ctx.bufferedFmState, ctx.scarPilotI, ctx.scarPilotQ,
            ctx.state_scar_pll, ctx.if_fs, ctx.mixedStereoFilt, ctx.expander, ctx.mono_audio_block, ctx.audio_coeff, ctx.state_audio_mono, ctx.audio_decim,
            ctx.phase_offset_audio, ctx.stereoAudioBlock, ctx.stereo_audio_coeff, ctx.state_stereo_base, ctx.phase_offset_stereo, ctx.block_id, ctx.audio_Fs, ctx.pcm_buf);
            
        auto end = std::chrono::high_resolution_clock::now();
        auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        benchmark::DoNotOptimize(duration_us);

        benchmark::ClobberMemory();
    }
}
BENCHMARK(Bench_Sequential)->RangeMultiplier(2)->Range(MIN_BLOCK_SIZE, MAX_BLOCK_SIZE);

// --- Test 4: True Pipelined Multithreading ---
static void Bench_Pipelined(benchmark::State& state) {
    SDR_Context ctx(state.range(0));
    ctx.block_id = 100;
    
    for (auto _ : state) {
        ctx.reset_states();

        auto start = std::chrono::high_resolution_clock::now();

        std::thread t_rf(rf_thread, ctx.block_size, ctx.block_id, std::ref(ctx.block_data), std::ref(ctx.i_b_half), 
            std::ref(ctx.q_b_half), std::ref(ctx.i_ds), ctx.rf_coeff, std::ref(ctx.state_i_lpf_100k),
            ctx.rf_decim, ctx.phase_offset_i, std::ref(ctx.q_ds), std::ref(ctx.state_q_lpf_100k), ctx.phase_offset_q,
            std::ref(ctx.fm_demod), std::ref(ctx.state_I_last_sample), std::ref(ctx.state_Q_last_sample));

        std::thread t_audio(audio_thread, ctx.format, std::ref(ctx.scarFilt), ctx.fm_demod_pipeline_in, ctx.scar_coeff, std::ref(ctx.state_scar_pilot), ctx.scar_phase_offset,
            std::ref(ctx.schFilt), ctx.sch_coeff, std::ref(ctx.state_sch_stereo), ctx.sch_phase_offset, std::ref(ctx.bufferedFm_demod), std::ref(ctx.bufferedFmState), std::ref(ctx.scarPilotI), std::ref(ctx.scarPilotQ),
            std::ref(ctx.state_scar_pll), ctx.if_fs, std::ref(ctx.mixedStereoFilt), ctx.expander, std::ref(ctx.mono_audio_block), ctx.audio_coeff, std::ref(ctx.state_audio_mono), ctx.audio_decim,
            ctx.phase_offset_audio, std::ref(ctx.stereoAudioBlock), ctx.stereo_audio_coeff, std::ref(ctx.state_stereo_base), ctx.phase_offset_stereo, ctx.block_id, ctx.audio_Fs, ctx.pcm_buf);

        t_rf.join();
        t_audio.join();

        auto end = std::chrono::high_resolution_clock::now();
        auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        benchmark::DoNotOptimize(duration_us);

        std::swap(ctx.fm_demod, ctx.fm_demod_pipeline_in);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(Bench_Pipelined)->RangeMultiplier(2)->Range(MIN_BLOCK_SIZE, MAX_BLOCK_SIZE);