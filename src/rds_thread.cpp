#include "../include/dy4.h"
#include "../include/rds_thread.h"
#include "../include/RDSFrameSynchro.h"
#include "../include/RDSApplicationLayer.h"
#include "../include/rds_helpers.h"

RDSBlockOffset convert_offset(Offset off) {
    switch (off) {
        case Offset::A: return A;
        case Offset::B: return B;
        case Offset::C:
        case Offset::Cp: return C;
        case Offset::D: return D;
        default: return A;
    }
}

void rds_thread(int mode, int block_id, float block_duration, const std::vector<real>& fm_demod, //INPUT
                std::vector<real>& rds_coeff, std::vector<real>& state_rds, std::vector<real>& bufferedFM_state, int& band_phase_offset, //Channel Extraction
                float if_FS, std::vector<real>& car_rec_coeff, std::vector<real>& state_carr, std::vector<real>& state_RDS_PLL, int& rec_carr_p_off, //CARRIER RECOVERY
                std::vector<real>& demod_coeff, std::vector<real>& rrc_coeff, std::vector<real>& state_mix, std::vector<real>& state_mixQ, //DEMODULATION
                std::vector<real>& state_RRC, std::vector<real>& state_RRC_Q, int& phase_offset_RRC, int& phase_offset_RRC_Q, int& mixed_phase_offs, int& mixed_phase_offs_Q, //DEMODULATION
                int rds_expander,  int rds_decimator, //Resampling
                float& prev_symbol, float& prev_sample, float& cdr_offset, int& next_init_pos, float& state_last_symbol, //CDR
                std::vector<real>& state_leftover_symbols, bool& unlocked, int& manchesterMode, //MANCHESTER DECODING
                int& state_last_bit, //DIFF DECODER
                std::vector<int>& state_leftover_decoded_bits, RDSFrameSynchro& rds_synchro, // FRAME SYNC
                RDSApplicationLayer& rds_app // APP LAYER
            )
{

    std::vector<real> rds_band; //created everytime the thread runs
    std::vector<real> bufferedRDS_band; //created everytime the thread runs, only used inside the thread
    std::vector<real> RDSPilotI, RDSPilotQ; //these are created inside the thread by carrier recovery, adn only used in demodulate.
    std::vector<real> rrcFilt, rrcFiltQ; //outputted of demod, only rrcFilt is used in CDR, FiltQ never used.
    std::vector<real> symbols; //this is our outut from CDR and input the MD, otherwise it is not needed, just need to carry over the leftvover.
    std::vector<int> bits, new_bits; //output from MD and input to DD, afeyr that it is not used. new_bits is the output of DD, which is imedialy combined into combined_bits, so dont need to be apssed in.

    int rdsStart_blk; 
    float lockDuration = 5e-1;

    rdsStart_blk = static_cast<int>(std::ceil(lockDuration/block_duration));

    int sps;
    if (mode  == 0)
    {
        sps = 16;
    }
    else if (mode == 2){
        sps = 31;
    }
    else{
         throw std::runtime_error("Invalid mode in rds_thread");
    }

    //TRACKED INPUTS //Channel   mode, fm_demod, rds_Coeff, state_rds, bufferedFM_State, band_phase_offset
                     //Carrier   if_FS, car_rec_coeff, state_car, state_RDS_PLL, rec_carr_p_off
                     //DEMOD     demod_coeff, rrc_coeff, state_mix, state_mixQ, state_RRC, state_RRC_Q, phase_offset_RRC, phase_offset_RRC_Q, mixed_phase_offs
                     //CDR       block_id, next_init_pos, symbols, sps, prev_sample, prev_symbol, cdr_offset
                     //MD        state_leftover_symbols, unlocked, manchesterMode, state_last_symbol, bits
                     //DD        state_last_bit, new_bits
                     //Sync      state_leftover_decoded_bits RDSFrameSynchro &rds_synchro
                     //APP       RDSApplicationLayer &rds_app

    // 1. Channel Extraction
    rds_band = optimizedFilt(fm_demod, rds_coeff, state_rds, 1, band_phase_offset); 
    

    // thread delay block and carrier recovery

    // 2. Carrier Recovery
    std::thread delayblock([&](){
                bufferedRDS_band = delayBlock(rds_band, bufferedFM_state); // to match delay of carrier recoevry
            });
    std::thread recovery([&](){
                rds_carrier_recovery(rds_band, if_FS, car_rec_coeff, state_carr, state_RDS_PLL, RDSPilotI, RDSPilotQ, rec_carr_p_off); //outputs RDSPilotI and RDSPilotQ
            });         
    delayblock.join();
    recovery.join();       
 
    // 3. Demodulation & Resampling
    rds_demodulate(bufferedRDS_band, RDSPilotI, RDSPilotQ, demod_coeff, rrc_coeff, rds_expander, rds_decimator, state_mix,  //RRC filter coeff passed in
                    state_mixQ, state_RRC, state_RRC_Q, phase_offset_RRC, phase_offset_RRC_Q, rrcFilt, rrcFiltQ, mixed_phase_offs, mixed_phase_offs_Q); //rrcFilt, rrcFiltQ outputted


    if (block_id == rdsStart_blk) { //set this to be roughly 0.5 seconds

        float locMax;
        int maxPos;

        std::vector<real> window(rrcFilt.begin(),rrcFilt.begin() + sps + sps/2); //reason we use sps + sps/2 is to guarantee our window includes a peak or a trough
        findLocMax(window, locMax, maxPos); 
        next_init_pos = maxPos;
    }

    if (block_id >= rdsStart_blk){
        float errorGain = 0.001; // This is the proportional gain for the Timing Error Detection- Error function. 
        // Adjust as necessary by trial and error

        clock_data_recovery(rrcFilt, symbols, sps, next_init_pos, prev_sample, prev_symbol, cdr_offset, errorGain); //input the rrcFilter, output the symbols

        manchesterDecode(symbols, state_leftover_symbols, unlocked, manchesterMode, state_last_symbol, bits); //input symbols and outputs bits

        diffDecode(bits,state_last_bit, new_bits); // we pass in the state_last_bit because we XOR with prev bit, and we pass in 
        
        std::vector<int> combined_bits;
        if ((int)state_leftover_decoded_bits.size() > 0){
            combined_bits = state_leftover_decoded_bits; //combining the new bits and what was left from previous 
            state_leftover_decoded_bits.clear(); //clear becauser we will reuse it after block processing
        }
        combined_bits.insert(combined_bits.end(), new_bits.begin(), new_bits.end());

        // RDS BLOCK PROCESSING

        auto [valid_blocks, remaining_bits] = rds_synchro.process_block(combined_bits); //RDSFrameSynchro &rds_synchro NEEDS TO BE PASSED INTO THREAD
            //broken into the two parts
        state_leftover_decoded_bits.insert(state_leftover_decoded_bits.begin(),remaining_bits.begin(), remaining_bits.end());  //update leftovers with remaining from call

        //RDS APPLCATION LAYER
        //each valid block contains 26 bits (the block), block type (A,B,C,D)
        for (auto &block_pair : valid_blocks) { //through all valid blocks, each blockpair is vector<int>, offset (letter)

            const std::vector<int> &block_bits = block_pair.first; //the firt part is the block bits
            Offset off = block_pair.second; //the second is the offset letter

            RDSBlockOffset app_off = convert_offset(off); //convert the enum, Offset:: A-> A

            rds_app.add_block(block_bits, app_off); //give RDS block and process it
        }
        if (!valid_blocks.empty()) { //display blocks if we have some valid ones
                rds_app.display();
        }
    }
}
