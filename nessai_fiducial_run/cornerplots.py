# Make corner plot from results in results_file_name, plotting only parameters in columns_to_plot (for ordering see all_labels)
#

import bilby
import json
import corner
import numpy as np
import matplotlib.pyplot as plt

results_file_name = "results_fiducial_nessai_zero-noise.json"

all_labels = ['$\\mathcal{M}$', '$q$', '$a_1$', '$a_2$', '$\\theta_1$', '$\\theta_2$', '$\\Delta\\phi$', '$\\phi_{JL}$', '$\\mathrm{DEC}$', '$\\mathrm{RA}$', '$\\theta_{JN}$', '$\\psi$', '$\\Lambda_1$', '$\\Lambda_2$', '$t_c$']
# luminosity_distance= 100.0 not in labels for some reason. Also note discrepancy between phi_12 in truths and delta_phi in labels
all_truths = [1.486, 0.9, 0.04, 0.01, 1.0264673717225983, 2.1701305583885513, 5.0962562029664955, 2.518241237045709, 0.2205292600865073, 3.952677097361719, 0.25, 2.6973435044499543, 1500, 750, -0.01]

columns_to_plot = [0,1,2,7,12]
truths = []
labels = []
for i in range(len(all_truths)):
	if i in columns_to_plot:
		truths.append(all_truths[i])
		labels.append(all_labels[i])


kwargs = dict(
            bins=50, smooth=0.9, label_kwargs=dict(fontsize=16),
            title_kwargs=dict(fontsize=16), color='#0072C1',
            truth_color='tab:orange', quantiles=[0.16, 0.84], truths=truths,
            levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
            plot_density=False, plot_datapoints=True, fill_contours=True,
            max_n_ticks=3, hist_kwargs=dict(density=True))

with open(results_file_name) as f:
	results = json.load(f)
print(results.keys())
print(results['parameter_labels'])
samples = np.asarray(results["samples"]["content"])

samples = samples[:,columns_to_plot]

print(len(samples[0]))
print(type(samples[0]))
print(len(samples))
print(labels)

figure = corner.corner(samples, labels = labels, **kwargs)



figure.savefig("nessai_corner_zero-noise")
