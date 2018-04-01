#---------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University        #
# Created 2017-02-01, Modified 2017-06-07                 #
# Generate application scripts for Brown's OSCAR system   #
#---------------------------------------------------------#
#!/usr/bin/env python

import numpy as np
from parameters import *  # # # 

# Parameters
queue = 'ccmb-condo'
cores = '4'
runtime = '12:00:00'  # 2-20h at 4 processes
memory = '32G'  # ~ 10-50G at 4 processes
memory_large = '32G'
comment = 'Gravel w/ mut+rec var'

# Descriptors
prefix = '1000GPP1'
pop_list = ParametersShared().pop_list  # # # 
chrm_list = range(1, 23)  # #
chrm_large = range(1, 6)  # # 

# Generate random global seeds
def msms_seed():
    return str(np.random.random_integers(np.iinfo(np.int32).max))

# Create individual batch job
with open('slurm_scripts/{0}_main.sh'.format(prefix), 'w') as s:
    # Write shell header
    s.write('#!/bin/sh\n')
    # Loop over simulations
    for pop in pop_list:
        for chrm in chrm_list:
            label = '{0}_M_{1}_{2}'.format(prefix, pop, chrm)
            # Write line to shell file
            s.write('sbatch ./slurm_scripts/{0}.txt;\n'.format(label))
            # Write batch script file
            with open('slurm_scripts/{0}.txt'.format(label), 'w') as f:
                f.write('#!/bin/bash'+'\n')
                f.write('#SBATCH --qos={0}'.format(queue)+'\n')
                f.write('#SBATCH -n {0}'.format(cores)+'\n')
                f.write('#SBATCH -t {0}'.format(runtime)+'\n')
                if chrm in chrm_large:
                    f.write('#SBATCH --mem={0}'.format(memory_large)+'\n')
                else:
                    f.write('#SBATCH --mem={0}'.format(memory)+'\n')
                f.write('#SBATCH -J {0}'.format(label)+'\n')
                f.write('time python main_apps_main.py {0} {1} {2}\n'.format(str(prefix), str(pop), str(chrm)))  # #
                f.write('# Notes: {0}'.format(comment)+'\n')
