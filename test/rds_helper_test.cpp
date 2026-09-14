#include "gtest/gtest.h"
#include "dy4.h"
#include <vector>
#include <cmath>
#include <algorithm>

#include "rds_helpers.h"

namespace {

    class RDS_Helper_Fixture : public ::testing::Test {
        public:
            const real EPSILON = 1e-4;

            std::vector<real> test_signal;
            std::vector<real> symbols;
            std::vector<int> bits;

            void SetUp() override {
                test_signal = {0.5, -0.6, 0.7, -0.8, 0.9, -1.0, 0.3, -0.2};
            }
    };

    TEST_F(RDS_Helper_Fixture, findLocMax_BASIC) { //simple test to find max of the data
        float locMax = 0;
        int maxPos = -1;

        findLocMax(test_signal, locMax, maxPos);

        EXPECT_NEAR(locMax, 1.0, EPSILON);
        EXPECT_EQ(maxPos, 5);
    }

    TEST_F(RDS_Helper_Fixture, CDR_AMPLITUDE_VARIATION) {

        std::vector<real> input = { //use values that arent perfect, how it wil be for our actual data 
            0.8,1.2,
            -1.1,-0.9,
            1.05,0.95,
            -0.7,-1.3
        };

        int sps = 2; //were gonna do 2 samples per symbol 
        int initPos = 0; // starting index in signal
        float prevSample = 0; // previous sample
        float prevSymbol = 0; //prev sybmol
        float offset = 0; // timing offset 

        std::vector<real> out;  // recovered symbolsample values

        clock_data_recovery(input, out, sps, initPos, prevSample, prevSymbol, offset);

        std::vector<int> expected = {1, -1, 1, -1}; // Expected symbol sequence (basic but good test)

        EXPECT_EQ(out.size(), expected.size());

        for (int i = 0; i < out.size(); i++) {
            int sign = (out[i] > 0) ? 1 : -1;
            EXPECT_EQ(sign, expected[i]);
        }
    }

    TEST_F(RDS_Helper_Fixture, diffDecode_INVERSION_CHECK) {
        std::vector<int> input = {1, 0, 1, 1};
        std::vector<int> inputOpp = {0, 1, 0, 0};
        std::vector<int> output;
        std::vector<int> outputOpp;
        int state_last_bit = 0;
        int state_last_bitOpp = 1;

        diffDecode(input, state_last_bit, output);
        diffDecode(inputOpp, state_last_bitOpp, outputOpp);

        EXPECT_EQ(output.size(), outputOpp.size()); //comapre output sizes
        EXPECT_EQ(output, outputOpp); //compare decoded outputs, should match since we apply 180 degree phase to Opp
        EXPECT_EQ(inputOpp.back(), state_last_bitOpp); //make sure last bit is updated correctly
        EXPECT_EQ(input.back(), state_last_bit);
    }

    TEST_F(RDS_Helper_Fixture, manchesterDecode_NOMRAL_SEQUENCE) {

        std::vector<int> expected_bits = {1,0,1,1,0,0,1,0,1,0};


        std::vector<float> input = { //mancehsyer encoding to -1 and 1 similuatijg the wave
            1,-1,   // 1
            -1,1,   // 0
            1,-1,   // 1
            1,-1,   // 1
            -1,1,   // 0
            -1,1,   // 0
            1,-1,   // 1
            -1,1,   // 0
            1,-1,   // 1
            -1,1    // 0
        };

        std::vector<float> leftover;
        std::vector<int> out_bits;

        bool unlocked = true;
        int mode = 0;
        float last_symbol = NAN;

        manchesterDecode(input, leftover, unlocked, mode, last_symbol, out_bits);

        EXPECT_EQ(out_bits.size(), expected_bits.size());

        EXPECT_EQ(out_bits, expected_bits);
    }
    TEST_F(RDS_Helper_Fixture, manchesterDecode_ODD_SEQUENCE) {

        std::vector<int> expected_bits = {0, 1,0,1,1,0,0,1,0,1,0}; //same black as last test, but starts with 0, 
                                                                   //since we have NAN replaced by 0, and the input starts with a 1                              
        std::vector<float> input = { //mancehsyer encoding to -1 and 1 similuatijg the wave
            1,
            1,-1,   // 1
            -1,1,   // 0
            1,-1,   // 1
            1,-1,   // 1
            -1,1,   // 0
            -1,1,   // 0
            1,-1,   // 1
            -1,1,   // 0
            1,-1,   // 1
            -1,1    // 0
        };

        std::vector<float> leftover;
        std::vector<int> out_bits;

        bool unlocked = true;
        int mode = 0;
        float last_symbol = -1; //since we are passing in a odd number of symbols we need to add the previous block last symbol

        manchesterDecode(input, leftover, unlocked, mode, last_symbol, out_bits);
        EXPECT_EQ(out_bits.size(), expected_bits.size());
        EXPECT_EQ(out_bits, expected_bits);
    }

    TEST_F(RDS_Helper_Fixture, rds_carrier_recovery_OUTPUT_CHECK) {

        int block_size = 100;
        std::vector<real> rds_band(block_size, 0.5);  // constant input
        float Fs = 240000.0;
        
        std::vector<real> car_rec_coeff = {0.1, 0.2, 0.4, 0.2, 0.1}; 
        
        std::vector<real> state_carr(car_rec_coeff.size() - 1, 0.0); 
        std::vector<real> state_RDS_pll = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0}; 

        std::vector<real> RDSPilotI;
        std::vector<real> RDSPilotQ;
        int rec_carr_p_off = 0;

        rds_carrier_recovery(
            rds_band, Fs, car_rec_coeff,
            state_carr, state_RDS_pll,
            RDSPilotI, RDSPilotQ,
            rec_carr_p_off
        );
        bool all_zero_I = true;
        bool all_zero_Q = true;

        for (int i = 0; i < block_size; i++) {
            if (std::abs(RDSPilotI[i]) > 1e-6) all_zero_I = false;
            if (std::abs(RDSPilotQ[i]) > 1e-6) all_zero_Q = false;
        }

        EXPECT_FALSE(all_zero_I);
        EXPECT_FALSE(all_zero_Q);
        // signal should not be constant
        real first_I = RDSPilotI[0];
        bool all_same_I = true;

        for (int i = 1; i < block_size; i++) {
            if (std::abs(RDSPilotI[i] - first_I) > 1e-6) {
                all_same_I = false;
                break;
            }
        }

        EXPECT_FALSE(all_same_I);

        // ================= STATE CHECK =================
        EXPECT_EQ(state_carr.size(), car_rec_coeff.size() - 1);
    }

    TEST_F(RDS_Helper_Fixture, rds_demodulate_DIMENSION_CHECK) {
        int block_size = 120;
        std::vector<real> bufferedRDS_band(block_size, 0.5);
        std::vector<real> RDSPilotI(block_size, 1.0);
        std::vector<real> RDSPilotQ(block_size, 0.0);
        
        std::vector<real> demod_coeff = {0.2, 0.6, 0.2}; // Size 3
        std::vector<real> rrc_coeff = {0.1, 0.2, 0.4, 0.2, 0.1}; // Size 5
        
        int expander = 2;
        int decimator = 3;
        
        std::vector<real> state_mix(demod_coeff.size() - 1, 0.0);
        std::vector<real> state_mixQ(demod_coeff.size() - 1, 0.0);
        

        int phase_taps = (rrc_coeff.size() + expander - 1) / expander;
        std::vector<real> state_RRC(phase_taps - 1, 0.0);
        std::vector<real> state_RRC_Q(phase_taps - 1, 0.0);
        
        int phase_offset_RRC = 0;
        int phase_offset_RRC_Q = 0;

        std::vector<real> rrcFilt;
        std::vector<real> rrcFiltQ;
        int mixed_phase_offs = 0;
        int  mixed_phase_offs_Q = 0;
        rds_demodulate(bufferedRDS_band, RDSPilotI, RDSPilotQ, 
                       demod_coeff, rrc_coeff, expander, decimator, 
                       state_mix, state_mixQ, state_RRC, state_RRC_Q, 
                       phase_offset_RRC, phase_offset_RRC_Q, 
                       rrcFilt, rrcFiltQ, mixed_phase_offs, mixed_phase_offs_Q);

        // Expected output size = ceil((block_size * expander) / decimator)
        // (120 * 2) / 3 = 80
        int expected_out_size = (block_size * expander - 0 /* initial phase_offset */ + decimator - 1) / decimator;
        
        EXPECT_EQ(rrcFilt.size(), expected_out_size);
        EXPECT_EQ(rrcFiltQ.size(), expected_out_size);

        EXPECT_EQ(state_mix.size(), demod_coeff.size() - 1);
        EXPECT_EQ(state_mixQ.size(), demod_coeff.size() - 1);

        EXPECT_EQ(state_RRC.size(), phase_taps - 1);
        EXPECT_EQ(state_RRC_Q.size(), phase_taps - 1);
    }

}