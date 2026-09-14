#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

def impulseResponseRootRaisedCosine(Fs, N_taps):

	"""
	Root raised cosine (RRC) filter

	Fs  		sampling rate at the output of the resampler in the RDS path
				sampling rate must be an integer multipler of 2375
				this integer multiple is the number of samples per symbol

	N_taps  	number of filter taps

	"""

	# duration for each symbol - should NOT be changed for RDS!
	T_symbol = 1/2375.0

	# roll-off factor (must be greater than 0 and smaller than 1)
	# for RDS a value in the range of 0.9 is a good trade-off between
	# the excess bandwidth and the size/duration of ripples in the time-domain
	beta = 0.90

	# the RRC inpulse response that will be computed in this function
	impulseResponseRRC = np.empty(N_taps)

	for k in range(N_taps):
		t = float((k-N_taps/2))/Fs
		# we ignore the 1/T_symbol scale factor
		if t == 0.0: impulseResponseRRC[k] = 1.0 + beta*((4/math.pi)-1)
		elif t == -T_symbol/(4*beta) or t == T_symbol/(4*beta):
			impulseResponseRRC[k] = (beta/np.sqrt(2))*(((1+2/math.pi)* \
					(math.sin(math.pi/(4*beta)))) + ((1-2/math.pi)*(math.cos(math.pi/(4*beta)))))
		else: impulseResponseRRC[k] = (math.sin(math.pi*t*(1-beta)/T_symbol) +  \
					4*beta*(t/T_symbol)*math.cos(math.pi*t*(1+beta)/T_symbol))/ \
					(math.pi*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

	# returns the RRC impulse response to be used by convolution
	return impulseResponseRRC


def polyResampler(xb, h, state, upS, downS, phase_offset):
	lenh = len(h)
	lenx = len(xb)

	phase_taps = (lenh + upS - 1)/ upS
	leny = (lenx * upS - phase_offset + downS - 1) / downS

	yb = np.zeros(int(leny))
	res = 0

	for i in range(int(leny)):
		n = phase_offset + i * downS
		phase = n % upS
		x_i = n //upS

		for k in range(int(phase_taps)):
			tap = phase + k*upS
			if tap >= lenh:
				break
			x_pos = x_i - k

			if x_pos >= 0 and x_pos < lenx:
				res+= h[int(tap)] * xb[int(x_pos)]

			else:
				state_idx = phase_taps - 1 + x_pos
				if state_idx >= 0 and state_idx < len(state):
					res += h[int(tap)] * state[int(state_idx)]
		yb[i] = res

	last_n = phase_offset + (leny - 1)*downS
	phase_offset = (last_n + downS) - (lenx * upS)
	if lenx >= phase_taps - 1:
		state = xb[int(lenx - (phase_taps - 1)): lenx].copy()
		

	return yb, state, phase_offset


if __name__ == "__main__":

	pass
