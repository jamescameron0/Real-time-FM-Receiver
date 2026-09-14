#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import fmPll

from fmSupportLib import fmDemodUnwrap, fmDemodArctan, fmPlotPSD
# for take-home add your functions

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5 
# add other settings for audio, like filter taps, ...
audio_taps = 101
audio_Fc = 16e3

scar_Fc_low = 18.5e3
scar_Fc_high = 19.5e3
sch_Fc_low = 22e3
sch_Fc_high = 54e3

block_duration_ms = 40e-3
# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs take-home (il_vs_th = 1)
il_vs_th = 0

def filter(Fs, Fc, N_taps):
	center = (N_taps-1)/2
	scale_factor = 0
	h = [0] * N_taps
	for k in range(N_taps):
		n = k - center
		h[k] = 0
		if n == 0:
			h[k] = 2*Fc/Fs
		else:
			h[k] = np.sin(2 * np.pi * Fc * n / Fs) / (np.pi * n)
		h[k] = h[k] * (1/2 - 1/2 * np.cos(2 * np.pi * k / (N_taps - 1)))
		scale_factor += h[k]
	for k in range(N_taps):
		h[k] /= scale_factor
	return h

def convolve(x,h):
	y = np.zeros(len(h) + len(x) - 1)
	for n in range(len(y)):
		for k in range(len(h)):
			if n - k >= 0 and n-k < len(x):
				y[n] += h[k] * x[n-k]
	return y


def blockconv(xb, h, state):
	lenh = len(h)
	yb = np.zeros(len(xb))
	for n in range(len(yb)):
		for k in range(lenh):
			if n-k >= 0:
				yb[n] += h[k] * xb[n - k]
			else:
				yb[n] += h[k] * state[n + lenh - 1 - k]
	new_state = xb[-(lenh-1):]
	return yb, new_state

def block_conv_processing(monoAudioData, h, block_size):
	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = monoAudioData.shape)

	# start at the first block (with relative position zero)
	position = 0

	# initial filter state - state is the size of the impulse response minus 1
	# we need two channels for the state (one for each audio channel)
	state = np.zeros(shape = (len(h) - 1, 2))

	while True:
		# filter both left and right channels
		filtered_data[position:position + block_size, 0], state[:, 0] = \
			blockconv(monoAudioData[position:position + block_size, 0], h, state[:,0])
		# the filter state has been saved only for the first channel above
		# you will need to adjust the code for the second channel below
		filtered_data[position:position + block_size, 1], state[:, 1] = \
			blockconv(monoAudioData[position:position + block_size, 1], h, state[:,1])

		position += block_size

		# the last incomplete block is ignored
		if position > (len(monoAudioData) - block_size):
			break

	return filtered_data


def fmdemodlab(I, Q, previous_I = 0.0, previous_Q = 0.0):
	
    fm_demod = np.empty(len(I))
 
    # iterate through each of the I and Q pairs
    for k in range(len(I)):
 
        # extract the current in-phase and quadrature samples
        current_I = I[k]
        current_Q = Q[k]
 
        iprime = current_I - previous_I
        qprime = current_Q - previous_Q
 
        # use the atan2 function (four quadrant version) to detect
        # the phase difference directly between consecutive samples
        fm_demod[k] = ((current_I * qprime - current_Q * iprime)/(np.square(current_I) + np.square(current_Q))) \
			if np.square(current_I) + np.square(current_Q) > 0 else 0
 
        # save the state of the current I and Q samples
        # to compute the next phase difference
        previous_I = current_I
        previous_Q = current_Q
 
    # return both the demodulated samples as well as the last I and Q values
    # (the last I and Q samples are needed to enable continuity for block processing)
    return fm_demod, previous_I, previous_Q

def delayBlock(input_block, state_block):
    output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
    state_block = input_block[-len(state_block):]

    return output_block, state_block

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/iq_samples/stereo_l0_r9.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	'''
	# IQ data is normalized between -1 and +1 in 64-bit double format
	iq_data = (np.float64(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 64-bit double format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")
	'''

	# coefficients for the front-end low-pass filter
	#rf_coeff = signal.firwin(audio_taps, audio_Fc / (rf_Fs/rf_decim/2), window=('hann'))
	rf_coeff = signal.firwin(rf_taps, rf_Fc / (rf_Fs / 2), window=('hann'))
	# coefficients for the filter to extract mono audio
	# to be updated by you during the in-lab session based on firwin
    # same principle as for audio_coeff (but different arguments, of course)
	audio_coeff = signal.firwin(audio_taps, audio_Fc / (rf_Fs/rf_decim/2), window=('hann'))

	scar_coeff = signal.firwin(audio_taps, [scar_Fc_low / (rf_Fs/rf_decim/2), scar_Fc_high / (rf_Fs/rf_decim/2)], pass_zero='bandpass', window=('hann'))
	sch_coeff = signal.firwin(audio_taps, [sch_Fc_low / (rf_Fs/rf_decim/2), sch_Fc_high / (rf_Fs/rf_decim/2)], pass_zero='bandpass', window=('hann'))
	stereo_audio_coeff = signal.firwin(audio_taps, 38e3 / (rf_Fs/rf_decim/2), window=('hann'))
		
	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = audio_Fs * block_duration_ms * rf_decim * audio_decim * 2
	block_size = int(block_size)
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps - 1)
	state_q_lpf_100k = np.zeros(rf_taps - 1)
	state_phase = 0
	state_I_last_sample = 0
	state_Q_last_sample = 0

	# states for audio filters for block-based processing
	state_audio_mono = np.zeros(audio_taps - 1)    # mono audio filter state
	state_scar_pilot = np.zeros(audio_taps - 1)    # pilot tone filter state
	state_sch_stereo = np.zeros(audio_taps - 1)    # stereo channel filter state
	state_stereo_base = np.zeros(audio_taps - 1)   # stereo baseband filter state
	# state for stereo carrier PLL: [integrator, phaseEst, feedbackI, feedbackQ, trigOffset]
	state_scar_pll = [0.0 ,0.0 ,1.0 ,0.0 ,1.0 ,0.0 ,0]  # initialized as 0, will become a list after first PLL call
	# state for buffered mono and stereo delay (both need same size for synchronization)
	bufferedFmState = np.zeros(int((audio_taps - 1) / 2))
	# audio buffer that stores all the audio blocks
	monoAudioData = np.array([]) # used to concatenate filtered blocks (audio data)
	# stereo audio buffers
	stereoLeft = np.array([])
	stereoRight = np.array([])

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count + 1) * block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# downsample the I/Q data from the FM channel
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
			iq_data[block_count * block_size:(block_count + 1) * block_size:2],
			zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
			iq_data[block_count * block_size + 1:(block_count + 1) * block_size:2],
			zi=state_q_lpf_100k)
		
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		# filter to extract the FM channel (I samples are even, Q samples are odd)
		fm_demod, state_I_last_sample, state_Q_last_sample = fmdemodlab(i_ds, q_ds, state_I_last_sample, state_Q_last_sample)
		scarFilt, state_scar_pilot = signal.lfilter(scar_coeff, 1.0, fm_demod, zi=state_scar_pilot)
		schFilt, state_sch_stereo = signal.lfilter(sch_coeff, 1.0, fm_demod, zi=state_sch_stereo)

		# extract the mono audio data through filtering
		# downsample audio data
		bufferedfm_demod, bufferedFmState = delayBlock(fm_demod, bufferedFmState)
		monoAudioFilt, state_audio_mono = signal.lfilter(audio_coeff, 1.0, bufferedfm_demod, zi=state_audio_mono)
		monoAudioBlock  = monoAudioFilt[::audio_decim]
		
		# stereo carrier PLL to recover the 19 kHz pilot tone
		# the PLL locks onto the stereo pilot at 19 kHz
		scarPilotI, scarPilotQ, state_scar_pll = fmPll.fmPll(scarFilt, 19e3, rf_Fs/rf_decim, state_scar_pll, ncoScale=2.0)
		
		# demodulate the stereo channel by mixing with the recovered 38 kHz carrier
		# the stereo channel (DSB around 38 kHz) is demodulated using the 2x pilot frequency
		mixedStereoFilt = 2 * scarPilotI[:-1] * schFilt
		# low-pass filter to extract the baseband stereo signal (L-R)
		stereoBaseFilt, state_stereo_base = signal.lfilter(stereo_audio_coeff, 1.0, mixedStereoFilt, zi=state_stereo_base)
		# downsample the stereo audio data
		stereoAudioBlock = stereoBaseFilt[::audio_decim]

		# concatenate mono and stereo blocks to build the stereo output
		stereoLeft = np.concatenate((stereoLeft, 0.5*(monoAudioBlock + stereoAudioBlock)))
		stereoRight = np.concatenate((stereoRight, 0.5*(monoAudioBlock - stereoAudioBlock)))

		print(stereoAudioBlock)
		print(monoAudioFilt)
		
		# Below commented lines used for plotting
		
		# if block_count > 10 and block_count < 15:  # show plots from block 10 onwards for debugging

		# 	# plot PSD of selected block after FM demodulation
		# 	# (for easier visualization purposes we divide Fs by 1e3 to imply the kHz units on the x-axis)
		# 	# (this scales the y axis of the PSD, but not the relative strength of different frequencies)
		# 	ax0.clear()
		# 	fmPlotPSD(ax0, stereoAudioBlock, (audio_Fs / audio_decim) / 1e3, subfig_height[0], \
		# 		'scarFilt')
		# 	# output binary file name (where samples are written from Python)
		# 	fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
		# 	# create binary file where each sample is a 32-bit float
		# 	fm_demod.astype('float32').tofile(fm_demod_fname)


		# 	# PSD for stereo left channel (current block)
		# 	fmPlotPSD(ax1, mixedStereoFilt, (audio_Fs / rf_decim) / 1e3, subfig_height[1], 'schFilt')
		# 	# PSD for stereo right channel (current block)
		# 	fmPlotPSD(ax2, stereoRight, (rf_Fs / rf_decim) / 1e3, subfig_height[2], 'scarPilot')
		
			
		# 	'''
		# 	# create binary file where each sample is a 64-bit double
		# 	fm_demod.astype('float64').tofile(fm_demod_fname)
		# 	'''
			
		# 	# save figure to file
		# 	fig.savefig("../data/fmStereoBlock" + str(block_count) + ".png")

		block_count += 1
		

	print('Finished processing all the blocks from the recorded I/Q samples')

	# combine stereo channels into a single array
	stereo_data = np.column_stack((stereoLeft, stereoRight))


	# First lines for debugging - Writes left and right individually
	# Write audio data to file

	#wavfile.write("../data/fmStereoLeft.wav", int(audio_Fs), stereoLeft.astype('float32'))
	#wavfile.write("../data/fmStereoRight.wav", int(audio_Fs), stereoRight.astype('float32'))	
	wavfile.write("../data/fmStereo.wav", int(audio_Fs), stereoRight.astype('float32'))

	# uncomment assuming you wish to show some plots
	# plt.show()