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

samples_file = 'samples_bright.json'

with open(samples_file) as f:
    posterior_samples = json.load(f)["posterior"]["content"]

luminosity_distance = np.array(posterior_samples['luminosity_distance'])
H0 = z * c / luminosity_distance
print(len(H0))

posterior_H0, bin_edges = np.histogram(H0, np.arange(H0_min, H0_max + H0_step, H0_step))

figure = plt.figure()
ax = figure.add_axes([0.15, 0.1, 0.8, 0.8])
plt.title('H0 posterior from 1 bright run')
ax.set_ylim([0, max(posterior_H0 * 1.1)])
plt.stairs(posterior_H0, bin_edges)

plt.show()
