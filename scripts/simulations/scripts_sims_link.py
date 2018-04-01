#--------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University       #
# Created 2017-02-01, Modified 2017-11-15                #
# Generate simulation scripts for Brown's OSCAR system   #
#--------------------------------------------------------#
#!/usr/bin/env python

import numpy as np

# Parameters
queue = 'ccmb-condo'
cores = '4'
runtime = '48:00:00'  # 3-7h at 2 processes
memory = '8G'  # ~ 3G at 2 processes
comment = 'GravelFinal2011'
comment = 'Jouganos2017-YCC'

# Descriptors
global_seed = 2017111503
job_list = range(0, 40)
np.random.seed(global_seed)

# Generate random global seed
def msms_seed():
    return str(np.random.random_integers(np.iinfo(np.int32).max))

# Create individual batch job
with open('slurm_scripts/{0}.sh'.format(global_seed), 'w') as s:
    # Write shell header
    s.write('#!/bin/sh\n')
    # Loop over simulations
    for job in job_list:
        # Generate identifiers
        label = '{0}_{1}'.format(global_seed, job)
        # Write line to shell file
        s.write('sbatch ./slurm_scripts/{0}.txt;\n'.format(label))
        # Write batch script file
        with open('slurm_scripts/{0}.txt'.format(label), 'w') as f:
            f.write('#!/bin/bash'+'\n')
            f.write('#SBATCH --qos={0}'.format(queue)+'\n')
            f.write('#SBATCH -n {0}'.format(cores)+'\n')
            f.write('#SBATCH -t {0}'.format(runtime)+'\n')
            f.write('#SBATCH --mem={0}'.format(memory)+'\n')
            f.write('#SBATCH -J {0}'.format(label)+'\n')
            f.write('time python main_sims_link.py {0} {1} {2}\n'.format(str(global_seed), str(job), str(msms_seed())))
            f.write('# Notes: {0}'.format(comment)+'\n')
