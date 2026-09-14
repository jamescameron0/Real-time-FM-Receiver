#include "RDSFrameSynchro.h"
#include <iostream>

std::vector<std::vector<int>> P = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
    {1, 0, 1, 1, 0, 1, 1, 1, 0, 0},
    {0, 1, 0, 1, 1, 0, 1, 1, 1, 0},
    {0, 0, 1, 0, 1, 1, 0, 1, 1, 1},
    {1, 0, 1, 0, 0, 0, 0, 1, 1, 1},
    {1, 1, 1, 0, 0, 1, 1, 1, 1, 1},
    {1, 1, 0, 0, 0, 1, 0, 0, 1, 1},
    {1, 1, 0, 1, 0, 1, 0, 1, 0, 1},
    {1, 1, 0, 1, 1, 1, 0, 1, 1, 0},
    {0, 1, 1, 0, 1, 1, 1, 0, 1, 1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
    {1, 1, 1, 1, 0, 1, 1, 1, 0, 0},
    {0, 1, 1, 1, 1, 0, 1, 1, 1, 0},
    {0, 0, 1, 1, 1, 1, 0, 1, 1, 1},
    {1, 0, 1, 0, 1, 0, 0, 1, 1, 1},
    {1, 1, 1, 0, 0, 0, 1, 1, 1, 1},
    {1, 1, 0, 0, 0, 1, 1, 0, 1, 1}
};
 

RDSFrameSynchro::RDSFrameSynchro() : state(Offset::SYNC) {}
RDSFrameSynchro::~RDSFrameSynchro() {}

RDSFrameSynchro::Syndrome RDSFrameSynchro::compute_syndrome(const std::vector<int>& window) {
    Syndrome syndrome;
    for (int j = 0; j < 10; ++j) {
        int xor_sum = 0;
        for (int i = 0; i < 26; ++i) {
            // Equiv to matrix mult, but dig logic to save compute
            xor_sum ^= (window[i] & P[i][j]);
        }
        syndrome[j] = xor_sum;
    }
    return syndrome;
}

Offset RDSFrameSynchro::match_syndrome(const Syndrome& syn) {
    for (auto const& [off, expected] : SYNDROMES) {
        if (syn == expected) return off;
    }
    return Offset::SYNC; 
}

std::pair<std::vector<std::pair<std::vector<int>, Offset>>, std::vector<int>> 
RDSFrameSynchro::process_block(std::vector<int>& bitstream) {
    std::vector<std::pair<std::vector<int>, Offset>> valid_blocks;
    size_t i = 0;

    while (i + 26 <= bitstream.size()) {
        std::vector<int> window(bitstream.begin() + i, bitstream.begin() + i + 26);
        Syndrome win_syn = compute_syndrome(window);
        Offset detected_off = match_syndrome(win_syn);

        if (state == Offset::SYNC) {
            if (detected_off != Offset::SYNC) {
                state = detected_off;
                valid_blocks.push_back({window, detected_off});
                i += 26;
            } else {
                i++;
            }
        } else {
            bool valid_next = false;
            if (state == Offset::B) {
                valid_next = (detected_off == Offset::C || detected_off == Offset::Cp);
            } else {
                if (state == Offset::A)  valid_next = (detected_off == Offset::B);
                if (state == Offset::C || state == Offset::Cp) valid_next = (detected_off == Offset::D);
                if (state == Offset::D)  valid_next = (detected_off == Offset::A);
            }

            if (valid_next) {
                valid_blocks.push_back({window, detected_off});
                state = detected_off;
                i += 26;
            } else {
                state = Offset::SYNC;
                i++;
            }
        }
    }
    
    std::vector<int> remaining(bitstream.begin() + i, bitstream.end());
    return {valid_blocks, remaining};
}
