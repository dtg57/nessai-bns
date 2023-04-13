# Take the results from a dark BNS population injection and find the posterior on H0.
#
# This assumes we know the EoS exactly. For each run, take posterior on lambda and convert to posterior on m_source using EoS.
# m_detector / m_source - 1 gives posterior on redshift
# H0 = c * z / d_L, so get posterior on H0 using posterior on z and posterior on d_L
# Multiply all posteriors on H0 together to get overall likelihood
#

import json
import math
import matplotlib.pyplot as plt
import numpy as np
import random
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

# The default prior used for lambdas is uniform between 0 and 5000
# This function converts that prior to a different uniform prior between lambda_min and lambda_max
# This is done by removing all samples from H0 for which lambda_1 or lambda_2 are outside of this range (does not normalise)
def correct_lambda_prior(lambda_min, lambda_max, lambda_1, lambda_2, H0):
    H0_corrected = []
    for i in range(len(H0)):
        if lambda_min < lambda_1[i] < lambda_max and lambda_min < lambda_2[i] < lambda_max:
            H0_corrected.append(H0[i])
    return np.array(H0_corrected)

# Folder containing JSON samples files (append / to this)
results_folder = '3G results/'
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
        len_before = len(H0_both)
        # correct H0_both to remove samples with lambdas outside of given range
        H0_both = correct_lambda_prior(0, 2000, lambda_1, lambda_2, H0_both)
        len_after = len(H0_both)
        print('H0 samples removed to correct lambda priors:', len_before - len_after)
        # add H0_both to the dictionary of samples
        H0_samples_all[samples_file] = H0_both

H0_min = -1000
H0_max = 1500
H0_step = 1
H0_n = int((H0_max - H0_min) / H0_step)

# plot some of the individual H0 posteriors
for run in H0_samples_all:
    if random.randint(0,4) == 4:
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

full_posterior_H0 = np.ones(H0_n)
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
    # normalise this. If sum is zero or NaN then skip entirely (all samples have been removed due to having too high lambda)
    if sum(posterior_H0) > 0:
        posterior_H0 = posterior_H0 / (sum(posterior_H0) * (H0_step))
        # add to full posterior [should this be a multiplication?]
        full_posterior_H0 = full_posterior_H0 * posterior_H0

# normalise
full_posterior_H0 = full_posterior_H0 / (sum(full_posterior_H0) * (H0_step))
print(full_posterior_H0)

figure = plt.figure('full')
ax = figure.add_axes([0.15, 0.1, 0.8, 0.8])
plt.title('full H0 posterior')
plt.xlabel('$H_0$ / $Km s^{-1} Mpc^{-1}$')
plt.ylabel('Probability density')
ax.set_ylim([0, max(full_posterior_H0) * 1.1])
plt.stairs(full_posterior_H0, bin_edges)
# plot the true value of H0 used by Bilby
plt.axvline(x=67.74, color='orange')
plt.savefig('full_H0_posterior_eosmethod_all_extendedrange', dpi=500)

plt.show()
