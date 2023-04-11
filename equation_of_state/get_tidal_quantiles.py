#
# This takes a run result and finds the median and one-sigma errors for 4 quantities: mass_1, mass_2, lambda_1, lambda_2.
# The method bilby.core.result.Result.get_one_dimensional_median_and_error_bar() is used to find these quantities.
# These are used for scatterplot of lambda vs mass, to compare to the EoS
#

import bilby
import json
from os import listdir
import re

# Runs are done with a uniform prior between 0 and 5000 on tidal parameters lambda.
# We reweight the posterior to a different uniform prior between lambda_min and lambda_max (these must be within 0 and 5000)
# This is achieved by removing all samples which have lambda_1 or lambda_2 outside of the [lambda_min, lambda_max] range
def reweight_lambda_prior(posterior, lambda_min, lambda_max):
    keys = posterior.keys()
    values = [[] for i in range(len(keys))]
    new_posterior = dict(zip(keys, values))
    for i in range(len(posterior['lambda_1'])):
        if lambda_min < posterior['lambda_1'][i] < lambda_max and lambda_min < posterior['lambda_2'][i] < lambda_max:
            for key in keys:
                new_posterior[key].append(posterior[key][i])
    return new_posterior


# append / to this
results_folder = "results/"
output_file_name = "lambda_mass_quantiles.txt"

all_files = listdir(results_folder)

merged_files = []

for json_file in all_files:
    if (json_file[-17:] == "merge_result.json"):
        merged_files.append(json_file)

output_csv_lines = [
    "run_number,lambda_1_median,lambda_1_plus,lambda_1_minus,lambda_2_median,lambda_2_plus,lambda_2_minus,mass_1_median,mass_1_plus,mass_1_minus,mass_2_median,mass_2_plus,mass_2_minus"]
output_csv_lines[0] += '\n'

for json_file in merged_files:
    f = open(results_folder + json_file)
    samples = json.load(f)
    posterior = samples["posterior"]["content"]
    len_before = len(posterior['lambda_1'])
    # reweight the posterior to a new, more restricted uniform prior in lambdas (by removing any samples with lambdas outside of range given)
    posterior = reweight_lambda_prior(posterior, 0, 2000)
    len_after = len(posterior['lambda_1'])
    # if more than 75% of samples have been removed then discard this run
    keep_run = (len_after / len_before) > 0.25
    print(json_file, len_before, len_before - len_after, keep_run)
    if keep_run:
        result = bilby.core.result.Result(posterior=posterior)
        lambda_1_quantiles = result.get_one_dimensional_median_and_error_bar('lambda_1')
        lambda_2_quantiles = result.get_one_dimensional_median_and_error_bar('lambda_2')
        # make sure these are source frame masses. The "mass_1" and "mass_2" are detector frame, so redshifted.
        mass_1_quantiles = result.get_one_dimensional_median_and_error_bar('mass_1_source')
        mass_2_quantiles = result.get_one_dimensional_median_and_error_bar('mass_2_source')
        number_arr = re.findall('(?<=_data)\d+(?=_0_analysis_)', json_file)
        if len(number_arr) != 1:
            print('error: run number not found')
        else:
            number = number_arr[0]
            csv_columns = [number, str(lambda_1_quantiles.median), str(lambda_1_quantiles.plus),
                           str(lambda_1_quantiles.minus), str(lambda_2_quantiles.median), str(lambda_2_quantiles.plus),
                           str(lambda_2_quantiles.minus), str(mass_1_quantiles.median), str(mass_1_quantiles.plus),
                           str(mass_1_quantiles.minus), str(mass_2_quantiles.median), str(mass_2_quantiles.plus),
                           str(mass_2_quantiles.minus)]
            output_csv_lines.append(','.join(csv_columns) + '\n')

f = open(output_file_name, "w")
f.writelines(output_csv_lines)
print(output_csv_lines)
f.close()
