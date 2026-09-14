#include "../include/dy4.h"
#include "../include/rf_thread.h"

void rf_thread(int block_size, int block_id, std::vector<real> &block_data, std::vector<real> &i_b_half, 
        std::vector<real> &q_b_half, std::vector<real> &i_ds, std::vector<real> rf_coeff, std::vector<real> &state_i_lpf_100k,
        int rf_decim, int phase_offset_i, std::vector<real> &q_ds, std::vector<real> &state_q_lpf_100k, int phase_offset_q,
        std::vector<real> &fm_demod, float &state_I_last_sample, float &state_Q_last_sample){

    readStdinBlockData(block_size, block_id, block_data);
    if ((std::cin.rdstate()) != 0)
    {
        std::cerr << "End of input stream reached" << std::endl;
        exit(1);
    }
    
    for (int i = 0, j = 0; i < (int)block_data.size(); i += 2, j++) {
        i_b_half[j] = block_data[i];
        q_b_half[j] = block_data[i+1];
    }

    // Filter and Decimate, faster due to reduced number if computations 
    std::thread i_filt_t([&](){
        i_ds = optimizedFilt(i_b_half, rf_coeff, state_i_lpf_100k, rf_decim, phase_offset_i);
    });
    std::thread q_filt_t([&](){
        q_ds = optimizedFilt(q_b_half, rf_coeff, state_q_lpf_100k, rf_decim, phase_offset_q);
    });

    i_filt_t.join();
    q_filt_t.join();

    //FM demodulator
    fmDemodArctan(fm_demod, i_ds, q_ds, state_I_last_sample, state_Q_last_sample);
}