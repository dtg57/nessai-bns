# Take the results from a dark BNS population injection and find the posterior on H0.
#
# This assumes we know the EoS exactly. For each run, take posterior on lambda and convert to posterior on m_source using EoS.
# m_detector / m_source gives posterior on redshift
# H0 = c * z / d_L, so get posterior on H0 using posterior on z and posterior on d_L
# Multiply all posteriors on H0 together to get overall likelihood
#

import json
import math
import matplotlib.pyplot as plt
import numpy as np
from os import listdir

# speed of light in km/s
c = 299792

EOS_fit_coeffs = {
    'MPA1': [0.000042879938, -0.001410127, 0.020350306, -0.144805729, 0.232370002, 2.359861131]
}

EOS = 'MPA1'


def tidal_to_mass(tidal, coeffs):
    ln_tidal = math.log(tidal)
    return coeffs[0] * ln_tidal ** 5 + coeffs[1] * ln_tidal ** 4 + coeffs[2] * ln_tidal ** 3 + coeffs[
        3] * ln_tidal ** 2 + coeffs[4] * ln_tidal + coeffs[5]


# Folder containing JSON samples files (append / to this)
results_folder = '3G/results/'
# JSON samples files to include
samples_files = listdir(results_folder)
posterior_samples_all = {}
H0_samples_all = {}

for samples_file in samples_files:
    print('analysing: ' + samples_file)
    with open(results_folder + samples_file) as f:
        posterior_samples = json.load(f)["posterior"]["content"]
        posterior_samples_all[samples_file] = posterior_samples
        luminosity_distance = np.array(posterior_samples['luminosity_distance'])
        # calculate H0 from the first star's data
        mass_1_detector = np.array(posterior_samples['mass_1'])
        lambda_1 = np.array(posterior_samples['lambda_1'])
        mass_1_source = [tidal_to_mass(tidal, EOS_fit_coeffs[EOS]) for tidal in lambda_1]
        redshift_1 = mass_1_detector / mass_1_source - 1
        H0_1 = redshift_1 * c / luminosity_distance
        # calculate H0 from the second star's data
        mass_2_detector = np.array(posterior_samples['mass_2'])
        lambda_2 = np.array(posterior_samples['lambda_2'])
        mass_2_source = [tidal_to_mass(tidal, EOS_fit_coeffs[EOS]) for tidal in lambda_2]
        redshift_2 = mass_2_detector / mass_2_source - 1
        H0_2 = redshift_2 * c / luminosity_distance
        # since the H0 likelihoods from mass_1 and mass_2 are not independent (came from the same sampling process), we must take the average before combining with other runs
        H0_both = (H0_1 + H0_2) / 2
        # add H0_both to the dictionary of samples
        H0_samples_all[samples_file] = H0_both

H0_min = -1000
H0_max = 1000
H0_step = 10
H0_n = int((H0_max - H0_min) / H0_step)

'''
# plot all of the individual H0 posteriors
for run in H0_samples_all:
	H0_samples = H0_samples_all[run]
	H0_samples_clipped = []
	n_below = 0
	n_above = 0
	for H0 in H0_samples:
		if H0 < H0_min:
			n_below += 1
		elif H0 > H0_max:
			n_above += 1
		else:
			H0_samples_clipped.append(H0)
	print(len(H0_samples), n_below, n_above)
	figure = plt.figure(run)
	plt.title(run)
	plt.hist(x=H0_samples_clipped, bins=100)
'''

full_posterior_H0 = np.zeros(H0_n)
for run in H0_samples_all:
    H0_samples = H0_samples_all[run]
    H0_samples_clipped = []
    n_below = 0
    n_above = 0
    for H0 in H0_samples:
        if H0 < H0_min:
            n_below += 1
        elif H0 > H0_max:
            n_above += 1
        else:
            H0_samples_clipped.append(H0)
    print(run, len(H0_samples), n_below, n_above)
    posterior_H0, bin_edges = np.histogram(H0_samples_clipped, np.arange(H0_min, H0_max + H0_step, H0_step))
    # normalise this
    posterior_H0 = posterior_H0 / (sum(posterior_H0) * (H0_step))
    # add to full posterior
    full_posterior_H0 = full_posterior_H0 + posterior_H0

figure = plt.figure('full')
ax = figure.add_axes([0.15, 0.1, 0.8, 0.8])
plt.title('full H0 posterior')
ax.set_ylim([0, max(full_posterior_H0 * 1.1)])
plt.stairs(full_posterior_H0, bin_edges)

plt.show()
