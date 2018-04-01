#-----------------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University                #
# Created 2015-07-23, Modified 2017-10-25                         #
# Methods for processing simulated or real-world genomic dataset  #
#-----------------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import numpy as np
import pandas as pd
from features import *

def filter_low_freq(data_shared, filter_low):
    """Filter low frequency variants"""
    # round singleton and doubletons to fixation
    # perform this for each population separately
    pos_low = dict()
    pos_high = dict()
    for p in range(len(data_shared['haplotypes'])):
        hap = data_shared['haplotypes'][p]
        n = hap.shape[0]
        DAC_list = np.sum(hap[:, :], axis=0)
        pos_low[p] = np.where(DAC_list <= filter_low)[0]
        # pos_high[p] = np.where(DAC_list >= n-filter_low)[0]
        pos_high[p] = np.where(DAC_list >= n)[0]  # unfolded
        for i in pos_low[p]:
            if DAC_list[i] > 0:
                data_shared['haplotypes'][p][:, i] = np.zeros(n, int)
        for i in pos_high[p]:
            if DAC_list[i] < n:
                data_shared['haplotypes'][p][:, i] = np.ones(n, int)
    # monomorphic variants across all populations
    mono_low = sorted(list(reduce(set.intersection, [set(pos_low[p]) for p in pos_low])))
    mono_high = sorted(list(reduce(set.intersection, [set(pos_high[p]) for p in pos_high])))
    mono_comb = sorted(mono_low + mono_high)
    data_shared['varid'] = np.delete(data_shared['varid'], mono_comb)
    data_shared['genpos'] = np.delete(data_shared['genpos'], mono_comb)
    data_shared['physpos'] = np.delete(data_shared['physpos'], mono_comb)
    data_shared['recrates'] = np.delete(data_shared['recrates'], mono_comb)
    for p in range(len(data_shared['haplotypes'])):
        data_shared['haplotypes'][p] = np.delete(data_shared['haplotypes'][p], mono_comb, 1)
    return data_shared

def get_data_shared_app(pos_file, rec_file, hap_files):
    # posfile
    pos_temp = pd.read_table(pos_file, sep=' ', header=None)
    pos_varid = np.array(pos_temp[1].as_matrix())
    pos_genpos = np.array(pos_temp[2].as_matrix())
    pos_physpos = np.array(pos_temp[3].as_matrix())
    pos_recrates = np.array(pd.read_table(rec_file, sep=' ', header=None)[0].as_matrix())
    # haplotypes
    hap_haplotypes = dict()
    for p, _ in enumerate(hap_files):
        hap_haplotypes[p] = pd.read_table(hap_files[p], sep=' ', header=None).as_matrix().transpose()
        # print np.shape(hap_haplotypes[p])
    # data_shared
    data_shared = dict()
    data_shared['varid'] = pos_varid
    data_shared['genpos'] = pos_genpos
    data_shared['physpos'] = pos_physpos
    data_shared['recrates'] = pos_recrates
    data_shared['haplotypes'] = hap_haplotypes
    return data_shared

def get_meta_shared_app(key='NA', label='NA', job='NA', pop='NA', chrm='NA', pos_file='NA', hap_files='NA'):
    meta_shared = dict()
    meta_shared['key'] = key
    meta_shared['label'] = label
    meta_shared['job'] = job
    meta_shared['chr'] = chrm
    meta_shared['pop'] = pop
    meta_shared['pos_file'] = pos_file
    meta_shared['hap_files'] = hap_files
    return meta_shared

def out_features_app(features, features_file):
    with open(features_file, 'w') as f:
        f.write('\t'.join(features_list + features_stat) + '\n')
        for i, d in enumerate(features.downsample_list):
            f.write('\t'.join([check_main(features.features_values[s][features.pop_i][i], 4, s) for s in features_list + features_stat]) + '\n')

def out_features_app_standardize(features, features_file):
    with open(features_file, 'w') as f:
        f.write('\t'.join(features_list + features_standardize) + '\n')
        for i, d in enumerate(features.downsample_list):
            f.write('\t'.join([check_main(features.features_values[s][features.pop_i][i], 4, s) for s in features_list + features_standardize]) + '\n')

def out_metadata_app(metadata, metadata_file):
    with open(metadata_file, 'w') as m:
        m.write('\t'.join(metadata_list) + '\n')
        m.write('\t'.join([check_main(metadata.metadata_values[s], 4, s) for s in metadata_list]) + '\n')

def get_data_shared_sim(simulations):
    data_shared = dict()
    data_shared['varid'] = simulations.varid_list
    data_shared['genpos'] = simulations.genpos_list
    data_shared['physpos'] = simulations.physpos_list
    data_shared['recrates'] = simulations.recrates_list
    data_shared['haplotypes'] = simulations.haplotypes_list
    return data_shared

def get_meta_shared_sim(key='NA', label='NA', job='NA', pop='NA', chrm='NA', pos_file='NA', hap_files='NA'):
    """Same function as meta_shared_app"""
    return get_meta_shared_app(key, label, job, pop, chrm, pos_file, hap_files)

def out_features_sim(features, f):
    for i, d in enumerate(features.downsample_list):
        f.write('\t'.join([check_main(features.features_values[s][features.pop_i][i], 4, s) for s in features_list + features_stat]) + '\n')

def out_features_sim_standardize(features, f):
    for i, d in enumerate(features.downsample_list):
        f.write('\t'.join([check_main(features.features_values[s][features.pop_i][i], 4, s) for s in features_list + features_standardize]) + '\n')

def out_metadata_sim(metadata, m):
    m.write('\t'.join([check_main(metadata.metadata_values[s], 4, s) for s in metadata_list]) + '\n')
