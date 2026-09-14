/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include "../include/dy4.h"
#include "../include/filter.h"


// ======================================================================================================================== //
// ============================== IMPULSE RESPONSE FUNCTIONS ============================================================== //
// ======================================================================================================================== //

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h)
{
	// Allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	int center = (int)((num_taps - 1)/2);
	double scale_factor = 0;
	int n;

	for (int k = 0; k < num_taps; k++){
		n = k - center;
		if (n == 0){
			h[k] = 2*Fc/Fs;
		} else {
			h[k] = (std::sin(2*PI*Fc*n/Fs))/(PI*n);
		}
		
		h[k] = h[k]*(1 - (std::cos(2*PI*k/(num_taps - 1))))/2; 
		scale_factor = scale_factor +  h[k];	
	}

	for (int k = 0; k < num_taps; k++){
		h[k] = h[k]/scale_factor;
	}
}

void impulseResponseBPF(real f_low, real f_high, real fs, int N_taps, std::vector<real> &h){
	h.clear();
	h.resize(N_taps, 0.0);

    int center = (N_taps - 1)/2;
    real f_mid = (f_low + f_high)/2;
    real scale_factor = 0;

    for (int k = 0; k < N_taps; k++){
        int n = k - center;
        if (n == 0)
            h[k] = 2*f_high/fs - 2*f_low/fs; //handle center tap
        else 
            h[k] = std::sin(2*PI*f_high*n/fs)/(PI * n) - std::sin(2*PI*f_low*n/fs)/(PI * n); //sinc based response
        h[k] = h[k] * (0.5 - 0.5*std::cos(2*PI*k/(N_taps - 1))); //Hann Window
        scale_factor = scale_factor + h[k] * std::cos(2*PI*n*f_mid/fs);
    }
    
    for (int k = 0; k < N_taps; k++)
    {
        h[k] /= scale_factor;
    }
}

void impulseResponseRootRaisedCosine(const real Fs, const int N_taps, std::vector<real>& impulseResponseRRC){

    double T_symbol = 1.0/2375.0;

    double beta = 0.9;

	impulseResponseRRC.resize(N_taps, 0.0);

    int k;
	for(k = 0; k<N_taps; k++){
        double t = (k-N_taps/2.0)/Fs;

        if(t == 0){
            impulseResponseRRC[k] = 1.0+beta*((4.0/PI)-1.0);
        }
        else if (t == -T_symbol/(4*beta) || t == T_symbol/(4*beta)) {
            impulseResponseRRC[k] = (beta/std::sqrt(2.0))*(((1.0+2.0/PI)*std::sin(PI/(4.0*beta)))+((1.0-2.0/PI)*std::cos(PI/(4.0*beta))));
        }
        else { 
            impulseResponseRRC[k] = (std::sin(PI*t*(1.0-beta)/T_symbol)+4.0*beta*(t/T_symbol)*std::cos(PI*t*(1.0+beta)/T_symbol))/(PI*t*(1.0-std::pow(4.0*beta*t/T_symbol,2.0))/T_symbol);
        }
    }
}

// ======================================================================================================================== //
// ============================== CONVOLUTION FUCNTIONS =================================================================== //
// ======================================================================================================================== //

std::vector<real> blockConvolve(const std::vector<real> &xb, const std::vector<real> &h, std::vector<real> &state){
	int lenh = h.size();
	std::vector<real> yb;
	yb.clear();
	yb.resize(xb.size(), 0.0);

	for (size_t n = 0; n < yb.size(); n++){
		for (int k = 0; k < lenh; k++){
			if (((n - k) >= 0 ) && ((n - k) < xb.size()))
			{
				yb[n] += (h[k]*xb[n-k]);} 
			else {
				yb[n] += h[k]*state[n + lenh - 1 -k ];}}}
	std::vector<real> new_state(xb.end() - (lenh - 1), xb.end());
	state = new_state;

	return yb;
}

void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) //used for testing comparisons
{
	// Allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size() + h.size() - 1, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (size_t n = 0; n < y.size(); n++){
		for (size_t k = 0; k < h.size(); k++){
			if ((n - k >= 0) && (n - k < x.size())){
				y[n] += (h[k] * x[n-k]);
			}
		}
	}
}

std::vector<real> optimizedFilt(const std::vector<real> &xb, const std::vector<real> &h,\
                                std::vector<real> &state, const int decim, int &phase_offset){ 
    const int lenh = (int)h.size();
    const int lenx = (int)xb.size();
    const int leny = (lenx - phase_offset + decim - 1) / decim;

    std::vector<real> yb(leny, 0.0f);

    // Only the first i_split samples have n < lenh and need state access.
    // After that, all taps land within xb — no branching needed.
    const int i_split = std::min(leny, (lenh - phase_offset + decim - 1) / decim);

    // Precompute unroll limit for the hot path (lenh is constant across all iterations)
    const int unroll_lim = lenh - (lenh % 4);

    // Part 1: samples that need state 
    for (int i = 0; i < i_split; i++) {
        const int n = phase_offset + i * decim;
        float res = 0.0f;

        // Precomputed boundary: k < k_xb uses xb, k >= k_xb uses state
        const int k_xb = std::min(n + 1, lenh);

        // State portion (only runs for small n, exits quickly)
        for (int k = k_xb; k < lenh; k++) {
            res += h[k] * state[lenh - 1 + n - k];
        }

        // xb portion with 4x unrolling
        const int xb_unroll_lim = k_xb - (k_xb % 4);
        int k = 0;
        for (; k < xb_unroll_lim; k += 4) {
            res += h[k]   * xb[n - k]
                 + h[k+1] * xb[n - k - 1]
                 + h[k+2] * xb[n - k - 2]
                 + h[k+3] * xb[n - k - 3];
        }
        for (; k < k_xb; k++) {
            res += h[k] * xb[n - k];
        }

        yb[i] = res;
    }

    // Part 2: pure xb samples
    for (int i = i_split; i < leny; i++) {
        const int n = phase_offset + i * decim;
        float res = 0.0f;

        int k = 0;
        for (; k < unroll_lim; k += 4) {
            res += h[k]   * xb[n - k]
                 + h[k+1] * xb[n - k - 1]
                 + h[k+2] * xb[n - k - 2]
                 + h[k+3] * xb[n - k - 3];
        }
        for (; k < lenh; k++) {
            res += h[k] * xb[n - k];
        }

        yb[i] = res;
    }

    const int last_n = phase_offset + (leny - 1) * decim;
    phase_offset = (last_n + decim) - lenx;

    // Simplified from clear + new_state + assign
    state.assign(xb.end() - (lenh - 1), xb.end());

    return yb;
}

void polyResampler(const std::vector<real> &xb, const std::vector<real> &h,\
                                std::vector<real> &state, const int upS, const int downS, int &phase_offset, std::vector<real> &yb, bool rds){
    //Input Frequency = upS * Fs
    // Output Frequency = upS/downS * Fs
    int lenh = h.size();
    int lenx = xb.size();
    //const int phase_taps = lenh/upS;
    const int phase_taps = (lenh + upS - 1) / upS;

    const int leny = (lenx * upS - phase_offset + downS - 1) / downS;

	yb.clear();
	yb.resize(leny, 0.0);        
    // The point is to avoid computations that we KNOW will be discarded at the output
    // Use polyphase techinque to do so.
    
    // n is the upsampled freuqnecy domain
    // i is the rationally resampled frequency domain (desired output, factor of upS/downS * Fs)
    // k shifts the sub-filter which is defined by  e_k[n] = h[k + upS*n]
    // Basically, the k will be incremented to simulate the zero padding 
    // the purpose of this is to select the "branch" of one sub filter, using tap variable and only do that one, saving compute
    for (int i = 0; i < leny; i++)
    {
        const int n = phase_offset + i * downS; //The index in the upsampled domain, which is to be decimated by downS afterwards
        const int phase = n % upS;
        const int x_i = n / upS; // The index of the original sampling dowmain

        real res = 0.0;

        for (int  k = 0; k < phase_taps; k++) // only need to compute for size of subfilter
        {
            const int tap = phase + k * upS; // index into h

            if (tap >= lenh) break;

            const int x_pos = x_i - k; // input sample index

            if (x_pos >= 0 && x_pos < lenx) {
                //Multiply and accumulate to compute the y output at rationally resampled frequency
                res += h[tap] * xb[x_pos];
            } else {
                // Out of range of current scope, get it from the previous state blcok
                const int state_idx = phase_taps - 1 + x_pos;
                if (state_idx >= 0 && state_idx < (int)state.size())
                    res += h[tap] * state[state_idx];
            }
        }
        yb[i] = (rds) ? res : res*upS; // i is the final rationally resampled domain
    }

    const int last_n  = phase_offset + (leny - 1) * downS;
    phase_offset  = (last_n + downS) - (lenx * upS);

    state.clear();
    //std::vector<real> new_state(xb.end() - (phase_taps - 1), xb.end());
    //state = new_state;
    if ((int)xb.size() >= phase_taps - 1) {
    state.assign(xb.end() - (phase_taps - 1), xb.end());
    } else {
        std::vector<real> combined(state);
        combined.insert(combined.end(), xb.begin(), xb.end());
        state.assign(combined.end() - (phase_taps - 1), combined.end());
    }
}

// ======================================================================================================================== //
// ============================== FM DEMOD FUCNTIONS ====================================================================== //
// ======================================================================================================================== //

void fmDemod(std::vector<real> &fm_demod, const std::vector<real> &I,\
            const std::vector<real> &Q, float &prev_I, float &prev_Q){
    float curr_I, curr_Q, i_p, q_p;

    fm_demod.resize(I.size());

    for (size_t k = 0; k < I.size(); k++)
    {
        curr_I = I[k];
        curr_Q = Q[k];

        i_p = curr_I - prev_I;
        q_p = curr_Q - prev_Q;

        fm_demod[k] = ((curr_I * q_p - curr_Q * i_p)/(std::pow(curr_I, 2) + std::pow(curr_Q, 2)));

        prev_I = curr_I;
        prev_Q = curr_Q;
    }
}

void fmDemodArctan(std::vector<real> &fm_demod, const std::vector<real> &I, \
            const std::vector<real> &Q, float &prev_I, float &prev_Q){

    float curr_I, curr_Q, real, imag;
    fm_demod.resize(I.size());
    
    for (size_t i = 0; i < I.size(); i++)
    {
        curr_I = I[i];
        curr_Q = Q[i];

        real = curr_I * prev_I + curr_Q * prev_Q;
        imag = curr_Q * prev_I - curr_I * prev_Q;

        fm_demod[i] = std::atan2(imag, real);
        
        prev_I = curr_I;
        prev_Q = curr_Q;
    }
}

// ======================================================================================================================== //
// ============================== FILTER INITIALIZATION FUNCTION ========================================================== //
// ======================================================================================================================== //

void initFilters(Filters& f, float rf_Fs, float rf_Fc, int rf_taps, float if_fs, float audio_Fc, int audio_taps, float expander,
                float stereo_Fc, float scar_Fc_low, float scar_Fc_high, float sch_Fc_low, float sch_Fc_high, int rds_expander) 
{
    impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, f.rf_coeff);  // RF front-end LPF

    impulseResponseLPF(if_fs * expander, audio_Fc, expander * audio_taps, f.audio_coeff);   // Mono audio LPF

    impulseResponseBPF(scar_Fc_low, scar_Fc_high, if_fs, audio_taps, f.scar_coeff); // Stereo filters
    impulseResponseBPF(sch_Fc_low, sch_Fc_high, if_fs, audio_taps, f.sch_coeff);

    impulseResponseLPF(if_fs * expander, stereo_Fc, expander * audio_taps, f.stereo_coeff); //stereo audio coeff

    impulseResponseBPF(54e3, 60e3, if_fs, audio_taps, f.rds_coeff); // RDS Channel extraction
    impulseResponseBPF(113.5e3, 114.5e3, if_fs, audio_taps, f.car_rec_coeff); //Carroer recovery
    impulseResponseLPF(if_fs, 3e3, audio_taps, f.demod_coeff); // Demod LPF
    impulseResponseRootRaisedCosine(if_fs * rds_expander, audio_taps * rds_expander, f.rrc_coeff); //RRC filter
}




