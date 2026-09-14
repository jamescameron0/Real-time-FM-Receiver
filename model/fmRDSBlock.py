#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Department of Electrical and Computer Engineering
# McMaster University
#
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import fmPll
import fmRRC

from fmSupportLib import fmDemodUnwrap, fmDemodArctan, fmPlotPSD
from enum import Enum, auto

class offset(Enum):
      SYNC = auto() # brute force
      
      # offset blocks
      A = auto()
      B = auto()
      C = auto()
      C_p = auto()
      D = auto()

class RDSFrameSynchro:
    def __init__(self):
        self.state = offset.SYNC  # Initialized to Start in brute force search

    def compute_syndrome(self, bits_window):
        #Multiply 26-bit vector by parity check matrix P ; Results in 1 x 10 offset code 
        bits = np.array(bits_window)
        code = np.empty((P.shape[0], P.shape[1]))
        syndrome = np.empty((P.shape[1]))
        for j in range(P.shape[1]): # 10
             # First do element wise & between the bitstream and each column of P
             code[:, j] = bits_window & P[:, j]
        for i in range(P.shape[1]):
            # Then do the xor reduction within each column vector, to squash the matrix to 1x10
            syndrome[i] =  np.logical_xor.reduce(code[:, i])

        return tuple(syndrome.astype(int))

    # Known syndromes for each offset type
    SYNDROMES = {
        offset.A:    (1,1,1,1,0,1,1,0,0,0),
        offset.B:    (1,1,1,1,0,1,0,1,0,0),
        offset.C:    (1,0,0,1,0,1,1,1,0,0),
        offset.C_p:  (1,1,1,1,0,0,1,1,0,0),
        offset.D:    (1,0,0,1,0,1,1,0,0,0),
    }

    # Correct block offset sequence. 
    NEXT_STATE = {
        offset.A:    offset.B,
        offset.B:    None, # B is None because it can either be C or C'
        offset.C:    offset.D,
        offset.C_p:  offset.D,
        offset.D:    offset.A,
    }

    def match_syndrome(self, syndrome):
        # off is offset.X... syn is the vector assoiated with it
        for off, syn in self.SYNDROMES.items():
            if np.array_equal(syndrome, syn):                
                # if match is found, return the offset the "state"
                return off
        return None
    
    
    # Main FSM
    def process_block(self, bitstream):
        
        # print(bitstream)

        valid_blocks = []
        i = 0
        #print(bitstream)
        while (i < len(bitstream)):
            window = bitstream[i: i + 26] 
            if len(window) < 26:  # not enough bits left for a full block
                return valid_blocks, bitstream[i:]
            win_syn = self.compute_syndrome(window)
            detected_off = self.match_syndrome(win_syn)

            if self.state == offset.SYNC: # if in brute search mode                
                if detected_off is not None:
                    # If syndrome is matched, the next state should be the next one after the detected offset
                    self.state = detected_off
                    print("Offset found at Bit Pos ", i, "Block Type ", self.state)
                    valid_blocks.append((window, detected_off))
                    i += 26 # Slide the window down a full block
                else: 
                    # If not, slide window down by one bit and repeat
                    i += 1

            else:  # if already synchronized
                if self.state == offset.B: 
                    # If state is B and the detected offset is either C or C'
                    # Confirm that it is a valid syndrom, and it can move on
                    valid_next = detected_off in (offset.C, offset.C_p)
                else:
                    expected_offset = self.NEXT_STATE[self.state] # Assign the expected offset based on current state
                    # For all other states, just have to check with the expected offset
                    valid_next = (detected_off == expected_offset)

                if valid_next:
                    # once confirmed that the syndrome is valid, append it to valid blocks
                    valid_blocks.append((window, detected_off))
                    self.state = detected_off
                    print("Valid next offset found at Bit Pos ", i, "Block Type ", self.state)
                    i += 26
                else:
                    self.state = offset.SYNC
                    print("Synchronization lost at bit position", i)
                    i += 1
    
        return valid_blocks, np.array([], dtype=int)
 

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5 

block_duration_ms = 40e-3

def bits_to_int(bits):
    result = 0
    for b in bits:
        result = (result << 1) | int(b)
    return result

PTY_TABLE = {
    0:  "None",          1:  "News",           2:  "Information",
    3:  "Sports",        4:  "Talk",           5:  "Rock",
    6:  "Classic Rock",  7:  "Adult Hits",     8:  "Soft Rock",
    9:  "Top 40",        10: "Country",        11: "Oldies",
    12: "Soft music",    13: "Nostalgia",      14: "Jazz",
    15: "Classical",     16: "Rhythm & Blues", 17: "Soft R&B",
    18: "Language",      19: "Religious Music",20:"Religious Talk",
    21: "Personality",   22: "Public",         23: "College",
    24: "Spanish Talk",  25: "Spanish Music",  26: "Hip Hop",
    29: "Weather",       30: "Emergency Test", 31: "Emergency",
}

class RDSApplicationLayer:

    def __init__(self):
        self.pi_code    = None
        self.pty        = None
        self.pty_name   = "N/A"
        self._ps_chars  = ['_'] * 8        # assembled one segment at a time
        self._group_buf = []               # accumulates (window, offset) for one group
        self.ps_name    = '________'

    def add_block(self, window, blk_offset):

        # If we receive an A block, start (or restart) a new group
        if blk_offset == offset.A:
            self._group_buf = [(window, blk_offset)]
            return

        # Discard the block if no group has been started yet
        if not self._group_buf:
            return

        self._group_buf.append((window, blk_offset))

        # A complete group has exactly 4 blocks ending with D
        if blk_offset == offset.D and len(self._group_buf) == 4:
            self._decode_group(self._group_buf)
            self._group_buf = []  # ready for the next group

    def display(self):
        import sys
        pi_str = f"{self.pi_code:04X}" if self.pi_code is not None else "N/A"
        print(f"  PI code: {pi_str}", file=sys.stderr)
        print(f"  Program type: {self.pty} – {self.pty_name}", file=sys.stderr)
        print(f"  Program service: {self.ps_name}", file=sys.stderr)

    def _info_bits(self, window):
        """Extract the 16 information bits from a 26-bit block window."""
        return window[:16]

    def _decode_group(self, group):

        # Extract 16-bit information words
        blk_a = self._info_bits(group[0][0])
        blk_b = self._info_bits(group[1][0])
        blk_d = self._info_bits(group[3][0])

        # Block A: PI code 
        self.pi_code = bits_to_int(blk_a)

        # Block B: group type, version, PTY, group-specific bits
        # Bit layout (MSB → LSB, indices 0-15):
        # [0:4] group type (4 bits)
        # [4] version (0 = A variant, 1 = B variant)
        # [6:11] PTY (5 bits)
        # [11:16] group-type-specific
        group_type = bits_to_int(blk_b[0:4])
        version    = int(blk_b[4])
        pty_val    = bits_to_int(blk_b[6:11])
        self.pty       = pty_val
        self.pty_name  = PTY_TABLE.get(pty_val, f"PTY {pty_val}")

        if group_type == 0 and version == 0:
            # Group 0A – Programme Service (PS) name
            # Block B bits[14:16]: 2-bit segment address (0-3 → pairs of chars)
            seg = bits_to_int(blk_b[14:16])          # 0, 1, 2, or 3
            # Block D: two PS characters
            c1 = bits_to_int(blk_d[0:8])
            c2 = bits_to_int(blk_d[8:16])
            self._ps_chars[seg * 2] = chr(c1) 
            self._ps_chars[seg * 2 + 1] = chr(c2) 
            self.ps_name = ''.join(self._ps_chars)
            print(f"[RDS] PS segment {seg}: "
                  f"'{self._ps_chars[seg*2]}{self._ps_chars[seg*2+1]}' → PS='{self.ps_name}'")

        else:
            # Other group types are received but not decoded
            print(f"[RDS] Group type {group_type}{'A' if version==0 else 'B'} received")


def clock_data_recovery(rrcWF, sps, init_pos, prev_sample, prev_symbol, offset_val, K = 0.001):
     #inputs:
     #rrcWF is the post RRC wave form
     #sps is the preset samples per symbol
     # K is the proportional gain
    # Mueller Muller algorithm for Timing Error Detection

     symbols = []
     sampled_pts = []

     offset = offset_val
     i = init_pos

     while i < len(rrcWF):
        curr_sample = rrcWF[i]
        curr_symbol = 1.0 if curr_sample > 0 else -1.0

        if prev_symbol != curr_symbol: 
            timing_error = prev_symbol * curr_sample - (curr_symbol * prev_sample) 
        else: 
            timing_error = 0.0  
            
        offset += timing_error * K
        offset = np.clip(offset, -sps//2, sps//2) 

        symbols.append(curr_sample)
        sampled_pts.append(i)

        prev_sample = curr_sample
        prev_symbol = curr_symbol

        step = int(round(sps + offset)) 
        step = max(1, min(step, 2 * sps))  
        i += step
     
     # Calculate where the index should start in the next block
     next_init_pos = i - len(rrcWF)

     return np.array(symbols), np.array(sampled_pts), prev_sample, prev_symbol, offset, next_init_pos


def findLocMax(sps, rrcWindow):
     # find the local maximum within this frame
    locMax = 0


    for i in range(len(rrcWindow)):
        pt = abs(rrcWindow[i])
        if pt > locMax:
             locMax = pt
             maxPos = i
    
    return locMax, maxPos


def delayBlock(input_block, state_block):
    output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
    state_block = input_block[-len(state_block):]

    return output_block, state_block


P = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0], \
             [0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], \
             [0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0], \
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 1, 1, 0, 1, 1, 1, 0, 0], [0, 1, 0, 1, 1, 0, 1, 1, 1, 0], \
             [0, 0, 1, 0, 1, 1, 0, 1, 1, 1], [1, 0, 1, 0, 0, 0, 0, 1, 1, 1], [1, 1, 1, 0, 0, 1, 1, 1, 1, 1], \
             [1, 1, 0, 0, 0, 1, 0, 0, 1, 1], [1, 1, 0, 1, 0, 1, 0, 1, 0, 1], [1, 1, 0, 1, 1, 1, 0, 1, 1, 0], \
             [0, 1, 1, 0, 1, 1, 1, 0, 1, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 1, 1, 1, 0, 1, 1, 1, 0, 0], \
             [0, 1, 1, 1, 1, 0, 1, 1, 1, 0], [0, 0, 1, 1, 1, 1, 0, 1, 1, 1], [1, 0, 1, 0, 1, 0, 0, 1, 1, 1], \
             [1, 1, 1, 0, 0, 0, 1, 1, 1, 1], [1, 1, 0, 0, 0, 1, 1, 0, 1, 1]])


if __name__ == "__main__":

    in_fname = "../data/iq_samples/samples2.raw"
    

    raw_data = np.fromfile(in_fname, dtype='uint8')
    iq_data = (np.float32(raw_data) - 128.0) / 128.0
    
    Fs = rf_Fs / rf_decim

    sps = 16
    expander = 19
    decimator = 120

    rf_coeff = signal.firwin(rf_taps, rf_Fc / (rf_Fs / 2), window=('hann'))

    block_size = audio_Fs * block_duration_ms * rf_decim * audio_decim * 2
    #block_size = 23999999
    block_size = int(block_size)
    block_count = 0

    state_i_lpf_100k = np.zeros(rf_taps - 1)
    state_q_lpf_100k = np.zeros(rf_taps - 1)
    state_I_last_sample = 0
    state_Q_last_sample = 0
    state_rds = np.zeros(rf_taps - 1)
    state_carr = np.zeros(rf_taps - 1)
    state_mix = np.zeros(rf_taps - 1)
    state_mixQ = np.zeros(rf_taps - 1)

    state_RDS_pll = [0.0 ,0.0 ,1.0 ,0.0 ,1.0 ,0.0 ,0] 
    bufferedFmState = np.zeros((rf_taps - 1)//2)
    state_RRC = np.zeros(rf_taps*expander - 1)
    state_RRC_Q = np.zeros(rf_taps* expander - 1)

    state_last_symbol = None
    state_last_bit = 0
    # STEP 1. BANDPASS 54kHz to 60kHz
    
    rds_coeff = signal.firwin(rf_taps, [54e3/(Fs/2), 60e3/(Fs/2)], pass_zero='bandpass', window=('hann')) #bandpass filter, 54kHz to 60kHz
    car_rec_coeff = signal.firwin(rf_taps, [113.5e3/(Fs/2), 114.5e3/(Fs/2)], pass_zero='bandpass', window=('hann')) #bandpass filter, 113.5kHz to 114.5kH
    demod_coeff = signal.firwin(rf_taps, 3e3 / (Fs / 2), window=('hann')) # Low Pass cutting off at 3 kHz
    rrc_coeff = fmRRC.impulseResponseRootRaisedCosine(Fs*expander, 101 * expander)


    rds_bitstream = np.array([], dtype=int)
    all_bits = np.array([])
    state_leftover_symbols = np.array([])

    phase_offset_RRC = 0
    phase_offset_RRC_Q = 0

    state_leftover_decoded_bits = np.array([], dtype=int)    
    subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
    plt.rc('figure', figsize=(7.5, 10))	# the size of the entire figure
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace = .6)

    rds_synchro = RDSFrameSynchro()
    rds_app_layer = RDSApplicationLayer()

    cdr_offset = 0.0
    next_init_pos = 0

    unlocked = True

    prev_sample = 0.0
    prev_symbol = 0.0

    while (block_count + 1) * block_size < len(iq_data):

        print("********Processing block", block_count, "********")

        # ================================================ RF Front end    =======================================================================
        i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff,1.0,iq_data[block_count*block_size:(block_count+1)*block_size:2], zi=state_i_lpf_100k)

        q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff,1.0,iq_data[block_count*block_size+1:(block_count+1)*block_size:2],zi=state_q_lpf_100k)

        i_ds = i_filt[::rf_decim] #240kHz
        q_ds = q_filt[::rf_decim]

        fm_demod, state_I_last_sample, state_Q_last_sample = fmDemodArctan(i_ds, q_ds, state_I_last_sample, state_Q_last_sample)

        # ========================================================================================================================================



        # ================================================ Channel Extraction  =======================================================================
        rds_band, state_rds = signal.lfilter(rds_coeff, 1.0, fm_demod, zi=state_rds)
        bufferedRDS_band, bufferedFmState = delayBlock(rds_band, bufferedFmState)

        # ========================================================================================================================================



        # ================================================ Carrier Recovery  =======================================================================

        # Squaring Non-Linearity: point-wise multiplication of each samlple with itself
        # The point is to square the input signal which will double the phase, resulting in no net phase shift 
        rds_band_squared = rds_band * rds_band 

        # BPF
        rds_recovered_carr, state_carr = signal.lfilter(car_rec_coeff, 1.0, rds_band_squared, zi=state_carr)

        # PLL with NCO to 57 kHz
        RDSPilotI, RDSPilotQ, state_RDS_pll = fmPll.fmPll(rds_recovered_carr, 114e3, Fs, state_RDS_pll, ncoScale=0.5, normBandwidth=0.001)
        # ========================================================================================================================================


        # ================================================ Demod  =======================================================================
        mixedRDSsignalI = RDSPilotI[:-1] * bufferedRDS_band 
        mixedFilt, state_mix = signal.lfilter(demod_coeff, 1.0, mixedRDSsignalI, zi=state_mix)

        mixedRDSsignalQ = RDSPilotQ[:-1] * bufferedRDS_band 
        mixedFiltQ, state_mixQ = signal.lfilter(demod_coeff, 1.0, mixedRDSsignalQ, zi=state_mixQ)

        # Rational REsampler - adjusts IF sample rate to a new rate that is an integer multiple of the symbol rate 
        # 16 SPS for Mode 0 and 31 SPS for Mode 2
        # Fs = 240 KSamples/sec, 
        # For mode 0
        # Then output sample rate = 16 * 2375 = 38000 Samples/sec
        # To determine expander and decimator between IF and output, 
        # gcd(240000, 38000) = 2000
        # Expansion Factor = 38000/2000 = 19, Decimation Factor = 240000/2000 = 120

        mixedFiltUps = np.zeros(len(mixedFilt) * expander)
        mixedFiltUpsQ = np.zeros(len(mixedFiltQ) * expander)


        # RRC - Goal is to reduce Inter-Symbol Interference (ISI)
        # Uses Nyquist Filter to produce an impulse response that is non zero at one specific sampling point and zero otherwise

        for i in range(len(mixedFilt)):
            mixedFiltUps[expander*i] = mixedFilt[i] 
        rrcFiltUps, state_RRC = signal.lfilter(rrc_coeff, 1.0, mixedFiltUps, zi=state_RRC)

        for i in range(len(mixedFiltQ)):
              mixedFiltUpsQ[expander*i] = mixedFiltQ[i] 
        rrcFiltUpsQ, state_RRC_Q = signal.lfilter(rrc_coeff, 1.0, mixedFiltUpsQ, zi=state_RRC_Q)

        rrcFilt = rrcFiltUps[::decimator]
        rrcFiltQ = rrcFiltUpsQ[::decimator]
        #===================================================================================================================================
        
        
        # Clock Data Recovery (CDR)
        if block_count == 13: # Hard Coded Based on calculation to start after roughly 0.5 seconds to allow PLL to lock
            locMax, next_init_pos = findLocMax(sps, rrcFilt[ :sps + sps // 2])


        if block_count >= 13: # let PLL lock for the first 0.5 seconds

            symbols = np.array([])
            Symbols = np.array([])
            bits = np.array([])

            sampled_idx = np.array([])
            # Must determine algorithm for symbol sampling
            # Cant use a fixed value, but must continuously monitor to combat drifts.
            errorGain = 0.01

            symbols, sampled_idx, prev_sample, prev_symbol, cdr_offset, next_init_pos = clock_data_recovery(
                rrcFilt, sps, next_init_pos, prev_sample, prev_symbol, cdr_offset, errorGain)           
                         
            valid_idx = sampled_idx[sampled_idx < len(rrcFiltQ)].astype(int)
            symbolsQ  = rrcFiltQ[valid_idx]
            symbols   = symbols[:len(valid_idx)]

            # # ================================================ Data Processing  =======================================================================
            # # at this point, the waveform is now translated into symbols (Symbol: H or L)
            # # Now this data must be converted from symbols to a bitstream, then these data frames must be identified to be structured into blocks
            # # Then these blocks are extracted, parsed, then the corresponding ascii characters are displayed

            # # Manchester and Differential Decoding : 
            # # Convert the stream of symbls recovered in previous step into a bitstream
            # # The data at transmission is Manchester Encoded, and the extracted bitstream has a rate of 1187.5 bits/sec (two symbols per bit)
            # # Differential Decoding ensures that phase inversion which can happen at the PLL-stage does not affect the interpretation of the data

            if block_count >= 0: 
                if state_last_symbol is not None:
                    Symbols = np.concatenate(([state_last_symbol], symbols)) # To make sure the boundary bit gets paired
                else:
                    Symbols = symbols
                
                # Must add Weak-High detection

                if len(state_leftover_symbols) > 0:
                    Symbols = np.concatenate((state_leftover_symbols, Symbols))
                    state_leftover_symbols = np.array([])
                
                pos = 0

                if block_count == 13 or unlocked: # 'Unlocked' means that there were too many 'desyncs' conssecutively
        
                    bits0 = np.array([])
                    bits1 = np.array([])

                    desync0 = 0
                    desync1 = 0
                    # # Manchester Decoding "check" phase

                    # Offset 0
                    k = 0
                    while k < len(Symbols) - 1 and k < 24:
                        if Symbols[k] > 0 and Symbols[k + 1] < 0: # HL
                            bits0 = np.append(bits0, 1)
                        elif Symbols[k] < 0 and Symbols[k + 1] > 0:
                            bits0 = np.append(bits0, 0)
                        else:
                            desync0 += 1
                        k += 2


                    # Offset 1
                    k = 1
                    while k < len(Symbols) - 1 and k < 23:
                        if Symbols[k] > 0 and Symbols[k + 1] < 0: # LH
                            bits1 = np.append(bits1, 1)
                        elif Symbols[k] < 0 and Symbols[k + 1] > 0:
                            bits1 = np.append(bits1, 0)
                        else:
                            desync1 += 1
                        k += 2
                        
                    # Set the phase accordingly and append the processed bits, then go processes normlly
                    manchesterMode = 0 if desync0 <= desync1 else 1
                    bits = bits0 if manchesterMode == 0 else bits1

                    pos =  24 - manchesterMode
                    unlocked = False

                    print("Manchester Mode:", manchesterMode)

                desync = 0
                
                
                while(pos < len(Symbols) - 1):
                    if Symbols[pos] > 0 and Symbols[pos + 1] < 0: # HL
                            bits = np.append(bits, 1)
                            pos += 2
                            desync = 0
                    elif Symbols[pos] < 0 and Symbols[pos + 1] > 0: # LH
                            bits = np.append(bits, 0)
                            pos += 2
                            desync = 0
                    else: 
                        desync += 1
                        pos += 2
                        if desync > 1:
                            unlocked = True
                            state_last_symbol = None         
                            # FIX: Prevent negative indexing from wiping the array
                            safe_idx = max(0, pos - 5)
                            state_leftover_symbols = Symbols[safe_idx:] 
                            print("Manchester Unlocked!")
                            break

                if pos == len(Symbols) - 1:
                    state_last_symbol = Symbols[-1]
                else:
                    state_last_symbol = None


                if len(bits) > 0:
                    bits = bits.astype(int)
                    new_bits = np.array([], dtype=int)

                    # # Differential Decoding
                    all_bits = np.concatenate(([state_last_bit], bits))
                    for i in range(1, len(all_bits)):
                        new_bits = np.append(new_bits, all_bits[i] ^ all_bits[i-1])



                    state_last_bit = bits[-1]
                    # state_last_bit, is not. 

                    if len(state_leftover_decoded_bits) > 0:
                        combined_bits = np.concatenate((state_leftover_decoded_bits, new_bits))
                        state_leftover_decoded_bits = np.array([], dtype=int)
                    else:
                        combined_bits = new_bits      

                    print("Combined Bits size: ", len(combined_bits))
                    valid_blks, state_leftover_decoded_bits = rds_synchro.process_block(combined_bits)
                
                    # RDS Application 
                    for blk_window, blk_offset in valid_blks:
                        rds_app_layer.add_block(blk_window, blk_offset)

                    # Display decoded RDS data once per processing block
                    rds_app_layer.display()
# ========================================================================================================================================
        block_count += 1

