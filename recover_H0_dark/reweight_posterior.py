import bilby
import matplotlib.pyplot as ppl
import numpy as np
import random
import math

# coefficients for quintic fits to EOSs
EOS_fit_coeffs = {
    'MPA1': [0.000042879938, -0.001410127, 0.020350306, -0.144805729, 0.232370002, 2.359861131]
}
# EOS to use
EOS = 'MPA1'

# convert tidal coefficient lambda to mass for a given equation of state (fit given in coeffs array)
def tidal_to_mass(tidal, coeffs):
    ln_tidal = math.log(tidal)
    return coeffs[0] * ln_tidal ** 5 + coeffs[1] * ln_tidal ** 4 + coeffs[2] * ln_tidal ** 3 + coeffs[
        3] * ln_tidal ** 2 + coeffs[4] * ln_tidal + coeffs[5]


# Initialise Bilby prior objects
lambda_prior = bilby.gw.prior.Uniform(name='lambda', minimum = 0, maximum = 2000)
M_prior = bilby.gw.prior.Uniform(name='chirp_mass', minimum = 1.43, maximum = 2.3)
q_prior = bilby.gw.prior.Uniform(name='mass_ratio', minimum = 0.6, maximum = 1)
dl_prior = bilby.gw.prior.UniformComovingVolume(name='luminosity_distance', minimum=20, maximum=200, unit='Mpc')
# speed of light in km/s
c = 299792
# limits on masses (these are additional constraints needed alongside the M_prior and q_prior)
m_min = 1.2
m_max = 2.3
# how many samples to take from these priors before binning 
n_samples = 10**6
# store H0 samples here
H0_samples = []

for i in range(n_samples):
	if i%10**4 == 0:
		print(i)
	dl_sample = dl_prior.sample()
	q_sample = q_prior.sample()
	M_sample = M_prior.sample()
	lambda_sample = lambda_prior.sample()
	# convert lambda to m_source
	m_source_sample = tidal_to_mass(lambda_sample, EOS_fit_coeffs[EOS])
	# the chirp_mass and mass_ratio priors actually encode for the two masses m_1 and m_2. Which of these masses we are actually sampling in is a matter of chance, equal probability
	if random.randint(0,1):
		# formula for m_1 in terms of q and chirp_mass
		m_detector_sample = (1+q_sample)**(1/5) * q_sample**(2/5) * M_sample
	else:
		# formula for m_2
		m_detector_sample = (1+q_sample)**(1/5) / q_sample**(3/5) * M_sample
	# if mass calculated is outside range then ignore sample
	if m_min < m_detector_sample < m_max:
		z_sample = m_detector_sample / m_source_sample - 1
		H0_sample = c*z_sample/dl_sample
		H0_samples.append(H0_sample)

ppl.hist(x = H0_samples, bins = 10000, range = (-1000,10000))
ppl.show()
