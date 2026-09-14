/*
   Comp Eng 3DY4 (Computer Systems Integration Project)

   Department of Electrical and Computer Engineering
   McMaster University
   Ontario, Canada
*/

// This file shows how to write convolution benchmark functions, using Google C++ benchmark framework.
// (it is based on https://github.com/google/benchmark/blob/main/docs/user_guide.md)

#include <benchmark/benchmark.h>
#include "utils_bench.h"
#include "dy4.h"
#include "iofunc.h"
#include "filter.h"

#define RANGE_MULTIPLIER 2
#define MIN_INPUT_SIZE 32768  // minimum signal length
#define MAX_INPUT_SIZE (4 * MIN_INPUT_SIZE)
#define MIN_FILTER_SIZE 101 // FIR filter length
#define MAX_FILTER_SIZE (1 * MIN_FILTER_SIZE)

const int lower_bound = -1; // Random signal generation bounds
const int upper_bound = 1;

const int UPS = 147;
const int DOWNS = 800; // Resampling factorss

// --- Reference: convolve then downsample ---
static void Bench_filtDecimate_reference(benchmark::State& state) {
    int N = state.range(0); // signal length
    int M = state.range(1); // filter length
    const int DECIM = 5;

    std::vector<real> x(N);
    std::vector<real> h(M);
    generate_random_values(x, lower_bound, upper_bound); // Generate random test signal and filter
    generate_random_values(h, lower_bound, upper_bound);

    for (auto _ : state) {
        //convolve
        std::vector<real> y_filt; //full conv
        convolveFIR(y_filt, x, h);

        // downsample
        std::vector<real> y_out; //downsample output
        for (int i = 0; i < (int)y_filt.size(); i += DECIM)
            y_out.push_back(y_filt[i]);

        benchmark::DoNotOptimize(y_out); //prevent optimization from complirer
    }
}

BENCHMARK(Bench_filtDecimate_reference)
    ->ArgsProduct({benchmark::CreateRange(MIN_INPUT_SIZE, MAX_INPUT_SIZE, RANGE_MULTIPLIER),
                   benchmark::CreateRange(MIN_FILTER_SIZE, MAX_FILTER_SIZE, RANGE_MULTIPLIER)});

// --- optimizedFilt ---
static void Bench_optimizedFilt(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);
    const int DECIM = 5;

    std::vector<real> x(N);
    std::vector<real> h(M);
    generate_random_values(x, lower_bound, upper_bound);
    generate_random_values(h, lower_bound, upper_bound);

    std::vector<real> bench_state(M - 1, 0.0); // State holds previous M-1 sample
    int phase_offset = 0;

    for (auto _ : state) {
        std::fill(bench_state.begin(), bench_state.end(), 0.0);  // Reset state each iteration 
        phase_offset = 0;

        std::vector<real> y = optimizedFilt(x, h, bench_state, DECIM, phase_offset);  // Optimizd filter+decim

        benchmark::DoNotOptimize(y); //prevent compiler optimizaitons
    }
}

BENCHMARK(Bench_optimizedFilt)
    ->ArgsProduct({benchmark::CreateRange(MIN_INPUT_SIZE, MAX_INPUT_SIZE, RANGE_MULTIPLIER),
                   benchmark::CreateRange(MIN_FILTER_SIZE, MAX_FILTER_SIZE, RANGE_MULTIPLIER)});

static void Bench_resample_reference(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);

    std::vector<real> x(N);
    std::vector<real> h(M);
    generate_random_values(x, lower_bound, upper_bound);
    generate_random_values(h, lower_bound, upper_bound);

    for (auto _ : state) {
        std::vector<real> x_up(N * UPS, 0.0); //upsample
        for (int i = 0; i < N; i++) x_up[i * UPS] = x[i];

        std::vector<real> y_filt; //convolve
        convolveFIR(y_filt, x_up, h);

        std::vector<real> y_out; //downsample
        for (int i = 0; i < (int)y_filt.size(); i += DOWNS)
            y_out.push_back(y_filt[i]);

        benchmark::DoNotOptimize(y_out);
    }
}

BENCHMARK(Bench_resample_reference)
    ->ArgsProduct({benchmark::CreateRange(MIN_INPUT_SIZE, MAX_INPUT_SIZE, RANGE_MULTIPLIER),
                   benchmark::CreateRange(MIN_FILTER_SIZE, MAX_FILTER_SIZE, RANGE_MULTIPLIER)});

static void Bench_polyResampler(benchmark::State& state) {
    int N = state.range(0);
    int M = state.range(1);

    std::vector<real> x(N);
    std::vector<real> h(M);

    bool rds = false;

    generate_random_values(x, lower_bound, upper_bound);
    generate_random_values(h, lower_bound, upper_bound);

    const int phase_taps = (M + UPS - 1) / UPS; // Number of taps per polyphase branch
    std::vector<real> bench_state(phase_taps - 1, 0.0); // holds the previuos samples for blocks
    int phase_offset = 0; //phase allignemnt across blocks
    std::vector<real> y;

    for (auto _ : state) {
        std::fill(bench_state.begin(), bench_state.end(), 0.0);
        phase_offset = 0;

        polyResampler(x, h, bench_state, UPS, DOWNS, phase_offset, y, rds); //polyphase resamplijng

        benchmark::DoNotOptimize(y);
    }
}

BENCHMARK(Bench_polyResampler)
    ->ArgsProduct({benchmark::CreateRange(MIN_INPUT_SIZE, MAX_INPUT_SIZE, RANGE_MULTIPLIER),
                   benchmark::CreateRange(MIN_FILTER_SIZE, MAX_FILTER_SIZE, RANGE_MULTIPLIER)});