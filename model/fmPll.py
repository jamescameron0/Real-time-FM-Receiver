#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

def fmPll(pllIn, freq, Fs, pll_states, ncoScale = 2.0, phaseAdjust = 0.0, normBandwidth = 0.01):

	"""

	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output

	phaseAdjust		float
					phase adjust to be added to the NCO output only

	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)

	state 			to be added

	"""

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp

	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	ncoOutI = np.empty(len(pllIn)+1)
	ncoOutQ = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = pll_states[0]
	phaseEst = pll_states[1]
	feedbackI = pll_states[2]
	feedbackQ = pll_states[3]
	ncoOutI[0] = pll_states[4]
	ncoOutQ[0] = pll_states[5]
	trigOffset = pll_states[6]

	# note: state saving will be needed for block processing
	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator

		# internal oscillator
		trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + phaseEst
		feedbackI = math.cos(trigArg)
		feedbackQ = math.sin(trigArg)
		ncoOutI[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
		ncoOutQ[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)

	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state

	pll_states = [integrator, phaseEst, feedbackI, feedbackQ, np.cos(trigArg*ncoScale), np.sin(trigArg*ncoScale), trigOffset]

	return ncoOutI, ncoOutQ, pll_states

	# for RDS add also the quadrature NCO component to the output

if __name__ == "__main__":

	pass
