#---------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University  #
# Created 2015-07-23, Modified 2017-06-07           #
# Main script for processing application data       # 
#---------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import sys
from parameters import *
from features import *
from main import *

if __name__ == '__main__':
    params_shared = ParametersShared()
    params_shared.downsample_sites = None
    label = str(sys.argv[1])
    pop = str(sys.argv[2])
    chrm = str(sys.argv[3])
    if len(sys.argv) == 5:
        params_shared.downsample_sites = float(sys.argv[4])
    job = 'M'  # place holder
    pos_file = '../../data/applications/1000GP_Phase1_post/{0}_chr{1}.pos'.format(label, chrm)
    rec_file = '../../data/applications/1000GP_Phase1_post/{0}_chr{1}_rec'.format(label, chrm)
    hap_files = dict()
    for p, q in enumerate(params_shared.pop_list):
        hap_files[p] = '../../data/applications/1000GP_Phase1_post/{0}_chr{1}_pop{2}.hap'.format(label, chrm, q)
    features_file = '../../data/applications/1000GP_Phase1_main/{0}_chr{1}_pop{2}_features.txt'.format(label, chrm, pop)
    metadata_file = '../../data/applications/1000GP_Phase1_main/{0}_chr{1}_pop{2}_metadata.txt'.format(label, chrm, pop)
    key = '{0}_chr{1}_pop{2}'.format(label, chrm, pop)
    meta_shared = get_meta_shared_app(key, label, job, pop, chrm, pos_file, hap_files)
    data_shared = filter_low_freq(get_data_shared_app(pos_file, rec_file, hap_files), params_shared.filter_low_freq_var)
    metadata = get_metadata(params_shared, meta_shared)
    features = get_features(params_shared, data_shared, meta_shared, hapbin=True)
    out_metadata_app(metadata, metadata_file)
    out_features_app(features, features_file)
