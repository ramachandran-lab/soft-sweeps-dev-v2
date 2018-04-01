#-------------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University            #
# Created 2015-07-23, Modified 2017-06-07                     #
# Main script for generating simulated standardizing dataset  #
#-------------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import sys
from parameters import *
from simulations import *
from features import *
from main import *

# time python main_standardize.py 2017021201 0 1234512345 1 1 1 

if __name__ == '__main__':
    params_sims = ParametersSimulations()
    params_sims.chromosome_length = params_sims.chromosome_length_stan
    params_sims.window_step_size_sim = params_sims.window_step_size_stan
    label = 'NA'
    job = 'NA'
    seed = 'NA'
    if len(sys.argv) >= 2:
        label = str(sys.argv[1])
        if len(sys.argv) >= 3:
            job = str(sys.argv[2])
            if len(sys.argv) >= 4:
                seed = str(sys.argv[3])
                np.random.seed(int(seed))
                if len(sys.argv) == 7:
                    params_sims.simulations_ne = int(sys.argv[4])
                    params_sims.simulations_hs = int(sys.argv[5])
                    params_sims.simulations_ss = int(sys.argv[6])
    for pop in params_sims.pop_list:
        features_file = '../../data/simulations/{0}/{1}_{2}_{3}_features.txt'.format(label, label, job, pop)
        metadata_file = '../../data/simulations/{0}/{1}_{2}_{3}_metadata.txt'.format(label, label, job, pop)
        with open(features_file, 'w') as f:
            with open(metadata_file, 'w') as m:
                f.write('\t'.join(features_list + features_standardize) + '\n')
                m.write('\t'.join(metadata_list) + '\n')
                i, j = 0, 0
                while j < params_sims.simulations_ne:
                    key = '{0}_{1}_{2}_{3}_{4}'.format(label, job, pop, 'NE', str(i))
                    simulations = NeutralModel(params_sims, key, pop, msms_seed()).print_msms().postprocess_msms()
                    if simulations.check:
                        meta_shared = get_meta_shared_sim(key, label, job, pop)
                        data_shared = filter_low_freq(get_data_shared_sim(simulations), params_sims.filter_low_freq_var)
                        out_metadata_sim(get_metadata(simulations, meta_shared), m)
                        out_features_sim_standardize(get_standardize(simulations, data_shared, meta_shared, hapbin=True), f)
                        j += 1
                    else:
                        metadata = get_metadata(simulations, get_meta_shared_sim(key, label, job, pop))
                        metadata.metadata_values['stype'] = 'neutral-dropped'
                        out_metadata_sim(metadata, m)
                    i += 1
