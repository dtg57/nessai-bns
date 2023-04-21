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
import random
from os import listdir

# speed of light in km/s
c = 299792

# coefficients for quintic fits to EOSs
EOS_fit_coeffs = {
    'MPA1': [0.000042879938, -0.001410127, 0.020350306, -0.144805729, 0.232370002, 2.359861131]
}
# EOS to use
EOS = 'MPA1'
# min and max lambda for reweighting the prior
lambda_min = 0
lambda_max = 2000


# convert tidal coefficient lambda to mass for a given equation of state (fit given in coeffs array)
def tidal_to_mass(tidal, coeffs):
    ln_tidal = math.log(tidal)
    return coeffs[0] * ln_tidal ** 5 + coeffs[1] * ln_tidal ** 4 + coeffs[2] * ln_tidal ** 3 + coeffs[
        3] * ln_tidal ** 2 + coeffs[4] * ln_tidal + coeffs[5]

# calculate the location of the given quantile (between 0 and 1) from the pdf.
# pdf is given as an array, so we also pass locations of bin_edges and the step
def get_quantile(pdf, quantile, bin_edges, step):
    cumsum = 0
    if (len(pdf) != len(bin_edges) - 1):
        print('error finding quantile')
        return False
    for i in range(len(bin_edges) - 1):
        if cumsum * step >= quantile:
            return bin_edges[i + 1]
        cumsum += pdf[i]

    return False

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

# boundaries to impose on H0 as well as step size when plotting histograms
H0_min = 0
H0_max = 1500
H0_step = 1
H0_n = int((H0_max - H0_min) / H0_step)
n_rejected = 0
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
    H0_both = []
    len_before = len(H0_1)
    for i in range(len(H0_1)):
        # reject any samples where lambdas are outside given range
        if lambda_min < lambda_1[i] < lambda_max and lambda_min < lambda_2[i] < lambda_max:
            # If both H0 are within bounds, then take an average and store. If only one is then keep that one only. If both are negative then reject sample
            H0_1_fine = H0_min < H0_1[i] < H0_max
            H0_2_fine = H0_min < H0_2[i] < H0_max
            if H0_1_fine:
                if H0_2_fine:
                    H0_both.append((H0_1[i] + H0_2[i]) / 2)
                else:
                    H0_both.append(H0_1[i])
            elif H0_2_fine:
                H0_both.append(H0_2[i])
    len_after = len(H0_both)
    print('H0 samples removed:', len_before - len_after, str((len_before - len_after) / len_before*100)+'%')
    # if less than 75% of samples have been removed due to reweighting priors and bounds imposed on H0, add H0_both to the dictionary of samples
    if (len_before - len_after) / len_before < 0.75:
        H0_samples_all[samples_file] = H0_both
    else:
        n_rejected += 1
        print('rejected:', samples_file)

print('total runs rejected:', n_rejected)

full_posterior_H0 = np.ones(H0_n)
n_used = 0
# runs to ignore because of poorly recovered parameters (often with very little support at H0 > 0)
bad_runs = []#['run34.json', 'run35.json', 'run45.json', 'run13.json', 'run18.json', 'run21.json']#["run2_08.json", "run1_36.json", "run1_13.json", "run2_18.json", "run2_19.json", "run2_36.json", "run2_39.json", "run2_40.json", "run2_41.json", "run2_29.json", "run2_17.json", "run2_08.json","run2_06.json","run2_12.json", "run1_12.json","run2_38.json","run2_46.json","run2_47.json","run2_13.json","run2_10.json"] #['run45.json', 'run18.json', 'run34.json', 'run35.json', 'run21.json', 'run13.json']

for run in H0_samples_all:
    H0_samples = H0_samples_all[run]
    # plot some of the individual H0 posteriors
    if random.randint(0,1) == 0:
        figure = plt.figure(run)
        plt.title(run)
        plt.hist(x=H0_samples, bins=100)
    posterior_H0, bin_edges = np.histogram(H0_samples, np.arange(H0_min, H0_max + H0_step, H0_step))
    # If this is a bad run skip. If sum is zero or NaN then skip, although this shouldn't be the case
    if sum(posterior_H0) > 0 and not run in bad_runs:
        # normalise
        posterior_H0 = posterior_H0 / (sum(posterior_H0) * (H0_step))
        # combine with full posterior by multiplying
        n_used += 1
        full_posterior_H0 = full_posterior_H0 * posterior_H0
        print(run, posterior_H0[1040:1090])
    else:
        print("rejected in second stage:", run)

print('runs used:', n_used)

# normalise
print(full_posterior_H0)
full_posterior_H0 = full_posterior_H0 / (sum(full_posterior_H0) * (H0_step))
print('total runs used:', n_used)
figure = plt.figure('full')
ax = figure.add_axes([0.15, 0.1, 0.8, 0.8])
plt.title('full H0 posterior')
plt.xlabel('$H_0$ / $Km s^{-1} Mpc^{-1}$')
plt.ylabel('Probability density')
ax.set_ylim([0, max(full_posterior_H0) * 1.1])
ax.set_xlim([-60,140])
plt.stairs(full_posterior_H0, bin_edges)
# plot the true value of H0 used by Bilby
plt.axvline(x=67.74, color='orange')
# calculate and plot the median and 1 sigma quantiles
median = get_quantile(full_posterior_H0, 0.5, bin_edges, H0_step)
quantile_lower = get_quantile(full_posterior_H0, 0.159, bin_edges, H0_step)
quantile_upper = get_quantile(full_posterior_H0, 0.841, bin_edges, H0_step)
print('median:', median)
print('lower quantile:', quantile_lower)
print('upper quantile:', quantile_upper)
plt.axvline(x=median, color='black')
plt.axvline(x=quantile_lower, color='black', ls='--')
plt.axvline(x=quantile_upper, color='black', ls='--')
plt.savefig('full_H0_posterior_3G_dark_with-quantiles_remove-below-zero', dpi=500)
'''
with open('out-dark.txt', 'w') as output:
    output.write(",".join(map(str, full_posterior_H0)))
'''
plt.show()
