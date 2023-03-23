# Plots lambda vs mass on a log-linear scale, with error bars in x and y, from data stored in csv_file_name
# plots all the NS EoSs stored in EOS_fit_params on the same graph
#

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math


# coefficients for the quintic fit for EoS, and the upper and lower mass limits for that fit (see Pacilio et al 2022 fig 1)
EOS_fit_params = {
    'MPA1' : {
        'coeffs' : [-3.1086, 25.356, -81.978, 131.558, -108.982, 45.168],
        'mass_max' : 2.3,
        'mass_min' : 1.2,
    },
    'SLY' : {
        'coeffs' : [-19.92705089, 145.3323963, -420.6281596, 603.7086089, -434.4840802, 133.7840435],
        'mass_max' : 2.04,
        'mass_min' : 1.2,
    },
    'GNH3' : {
        'coeffs' : [-30.80439489, 220.0863968, -624.8954323, 880.5481614, -620.2152793, 184.1083202],
        'mass_max' : 1.95,
        'mass_min' : 1.2,
    }
}

def mass_to_tidal_param(m, coeffs):
    # coeffs are for a exp(quintic) between mass and lambda
    return math.exp(coeffs[0]*m**5 + coeffs[1]*m**4 + coeffs[2]*m**3 + coeffs[3]*m**2 + coeffs[4]*m + coeffs[5])

def reduced_chi_square(mass_median, lambda_median, lambda_plus, lambda_minus, EOS):
    chi_square = 0
    n_points = 0
    for i in range(len(mass_median)):
        mass_data = mass_median[i]
        if mass_data < EOS['mass_max'] and mass_data > EOS['mass_min']:
            lambda_data = lambda_median[i]
            lambda_model = mass_to_tidal_param(mass_data, EOS['coeffs'])
            lambda_sigma = lambda_minus[i] if lambda_data > lambda_model else lambda_plus[i]
            chi_square += ((lambda_data - lambda_model) / lambda_sigma)**2
            n_points += 1
    return chi_square / n_points


csv_file_name = "Lambda_mass_quantiles_correct.csv"

csv_content = pd.read_csv(csv_file_name)
print(np.concatenate((csv_content.lambda_1_median.to_numpy(), csv_content.lambda_2_median.to_numpy())))
lambda_median = np.concatenate((csv_content.lambda_1_median.to_numpy(), csv_content.lambda_2_median.to_numpy()))
lambda_plus = np.concatenate((csv_content.lambda_1_plus.to_numpy(), csv_content.lambda_2_plus.to_numpy()))
lambda_minus = np.concatenate((csv_content.lambda_1_minus.to_numpy(), csv_content.lambda_2_minus.to_numpy()))

mass_median = np.concatenate((csv_content.mass_1_median.to_numpy(), csv_content.mass_2_median.to_numpy()))
mass_plus = np.concatenate((csv_content.mass_1_plus.to_numpy(), csv_content.mass_2_plus.to_numpy()))
mass_minus = np.concatenate((csv_content.mass_1_minus.to_numpy(), csv_content.mass_2_minus.to_numpy()))
print(len(mass_median))
print(len(lambda_median))
plt.errorbar(lambda_median, mass_median, xerr=[lambda_minus, lambda_plus], yerr=[mass_minus,mass_plus], fmt='o', elinewidth=1, ms=3)


print(mass_to_tidal_param(1.4, EOS_fit_params['GNH3']['coeffs']))
for EOS in EOS_fit_params:
    mass_EOS = np.arange(EOS_fit_params[EOS]['mass_min'], EOS_fit_params[EOS]['mass_max'], 0.005)
    lambda_EOS = [mass_to_tidal_param(m, EOS_fit_params[EOS]['coeffs']) for m in mass_EOS]
    plt.plot(lambda_EOS, mass_EOS, label=EOS)
plt.xlabel('$\Lambda$')
plt.ylabel('$M / M_{\odot}$')
plt.legend()
#
plt.xscale('log')
#plt.savefig('lambda_vs_mass_MPA1_SLY_GNH3_corrected', dpi=500)

for EOS in EOS_fit_params:
    print(EOS, reduced_chi_square(mass_median, lambda_median, lambda_plus, lambda_minus, EOS_fit_params[EOS]))

plt.show()