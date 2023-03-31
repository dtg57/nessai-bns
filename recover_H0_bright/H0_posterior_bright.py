# Calculate and plot posterior on H0 from a single bright run.
# We know the redshift of EM counterpart to infinite precision.
# Take the luminosity distance posterior from this run and do H0 = c * z / d_L to find posterior on H0.
#

import json
import math
import matplotlib.pyplot as plt
import numpy as np
from os import listdir

# speed of light in km/s
c = 299792
# redshift of EM counterpart
z = 0.02222

H0_min = 30
H0_max = 100
H0_step = 0.5
H0_n = int((H0_max - H0_min) / H0_step)

# folder where runs are stored, append / to this
samples_folder = 'results/'
samples_files = listdir(samples_folder)
posterior_samples = {}

for samples_file in samples_files:
    with open(samples_folder + samples_file) as f:
        posterior_samples[samples_file] = json.load(f)["posterior"]["content"]

H0_likelihood_all = {}
H0_likelihood_combined = np.ones(H0_n)
for samples_file in posterior_samples:    
    luminosity_distance = np.array(posterior_samples[samples_file]['luminosity_distance'])
    H0 = z * c / luminosity_distance
    print(len(H0))
    H0_likelihood, bin_edges = np.histogram(H0, np.arange(H0_min, H0_max + H0_step, H0_step))
    # normalise so likelihood integrates to 1
    H0_likelihood = H0_likelihood / (H0_step * sum(H0_likelihood))
    H0_likelihood_all[samples_file] = H0_likelihood
    H0_likelihood_combined = H0_likelihood_combined * H0_likelihood
    figure = plt.figure(samples_file)
    ax = figure.add_axes([0.15, 0.1, 0.8, 0.8])
    plt.title('H0 posterior from 1 bright run: ' + samples_file)
    ax.set_ylim([0, max(H0_posterior * 1.1)])
    plt.stairs(H0_posterior, bin_edges)

figure = plt.figure('all')
ax = figure.add_axes([0.15, 0.1, 0.8, 0.8])
plt.title('H0 posterior from 1 bright run: combined')
ax.set_ylim([0, max(H0_likelihood_combined * 1.1)])
plt.stairs(H0_likelihood_combined, bin_edges)

plt.show()
