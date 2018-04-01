#-------------------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University                  #
# Created 2015-07-23, Modified 2017-10-25                           #
# Calculates summary statistics from simulated or application data  #
#-------------------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import os, numpy as np
from conversions import *
import multiprocessing

from scipy import signal
import matplotlib.pyplot as plt
from scipy.stats import threshold
from sklearn.preprocessing import normalize

metadata_list = [
    'key',  # identify run
    # application
    'chr',  # identify chromosome
    'pop',  # identify population
    # simulation
    'stype',  # type of run
    'model',  # description
    'seed',  # seed of run
    # pre-simulation
    'mutrate',  # mutation rate
    'recrate',  # recombination rate
    'spop',  # dummy variable is 1 if sweep
    'stime',  # sweep initial time
    'scoeff',  # sweep strength coefficient
    'sfreq',  # sweep initial freq
    # post-simulation
    'complete',  # is sweep complete or incomplete
    'sweep_freq_bin',  # binned population final sweep frequency
    'sweep_freq_pop',  # population final frequency at sweep
    'sweep_fix_time',  # fixation time in past of sweep
    'sweep_dur_time',  # duration time in past of sweep
    'origins',  # number of original beneficial lineages
]

features_list = [
    # shared site metadaata
    'key',  # identify key
    'physpos',  # physical position of site
    # simulation site metadata
    'dist_physpos',  # physical distance of site to sweep
    'site_class',  # selective class of site
]

features_stat = [
    # window average of recrate
    # # 'R_val',  # recombination rate at site
    'segsites_win',  # number of segregating sites
    'recrates_win',  # window average of recrate at site
    # average nucleotide diversity (window-based)
    'theta_Pi_win',  # window average nucleotide diversity
    # # 'theta_S_win',  # window average nucleotide diversity
    'K_win',  # window haplotype diversity (Depaulis and Veuille 1998)
    # allele frequency statistics (SNP-based, cross-population)
    # # 'DAF_snp',  # derived allele frequency at site
    # 'DDAF_win',  # delta derived allele frequency, max of abs DDAF (Hofer et al. 2009, Grossman et al. 2010)
    # 'DDAF_mean',  # delta derived allele frequency, max of abs DDAF (Hofer et al. 2009, Grossman et al. 2010)
    # fixation index statistics (window-based, cross-population)
    'FST_win',  # mean of two comparisons, ratio-of-average Hudson's FST (1992), following Bhatia et al. (2103)
    # 'PBS_win',  # population branch statistic (Yi et al. 2010), based on above FST_win
    # site frequency spectrum statistics (window-based)
    'Taj_D_win',  # Tajima's D (1983), following Achaz et al. (2009)
    'FayWu_H_win',  # Fay and Wu's H (2000), following Achaz et al. (2009)
    'Zeng_E_win',  # Zeng et al.'s (2006), following Achaz et al. (2009)
    # haplotype frequency statistics (window-based)
    # # 'H1_win',  # average haplotype homozygosity (Garud et al. 2015)
    'H2_win',  # excluding most frequent haplotype (Garud et al. 2015)
    'H12_win',  # combining two most frequent haplotypes (Garud et al. 2015)
    'H2H1_win',  # ratio of above haplotype homozygosities (Garud et al. 2015)
    # linkage disequilibrium statistics (window-based)
    'ZA_win',  # (Rozas et al. 2001), average adjacent r^2 LD
    'Wall_B_win',  # (Wall 1999), based on congruent adjacent pairs
    # # 'Wall_Q_win',  # (Wall 1999), based on congruent adjacent pairs
    # extended haplotype homozygosity statistics (SNP-based)
    'iHS_win',  # iHS (unnormalized) (Voight et al. 2006, Grossman et al. 2010)
    # 'iHS_max',  # iHS (unnormalized) (Voight et al. 2006, Grossman et al. 2010)
    # 'D_iHH_win',  # delta iHH (unnormalized) (Voight et al. 2006, Grossman et al. 2010)
    # 'D_iHH_max',  # delta iHH (unnormalized) (Voight et al. 2006, Grossman et al. 2010)
    # cross-population extended haplotype homozygosity statistics (SNP-based)
    'XPEHH_win',  # XPEHH, max over two comparisons, (Sabeti et al. 2007, Grossman et al. 2010)
    # 'XPEHH_max',  # XPEHH, max over two comparisons, (Sabeti et al. 2007, Grossman et al. 2010)
    # delta site frequency spectrum estimators, based on theta_Pi (window-based, cross-population)
    # 'D_theta_S_win',  # delta theta_Pi, max over ratio, based on Innan and Kim (2008)
    'D_theta_Pi_win',  # delta theta_Pi, max over ratio, based on Innan and Kim (2008)
    # derived intra-allelic nucleotide diversity (SNP-based, cross-population)
    # # 'DIND_Pi_min',  # derived inter-alleic nucleotide diversity (Barreiro et al. 2009)
]

features_standardize = [
    'DAF_snp',  # derived allele frequency at site
    # 'DDAF_snp',  # delta derived allele frequency at site
    'iHS_snp',  # iHS (unnormalized) (Voight et al. 2006, Grossman et al. 2010)
    # 'D_iHH_snp',  # delta iHH (unnormalized) (Voight et al. 2006, Grossman et al. 2010)
    # 'DIND_Pi_snp',  # derived inter-alleic nucleotide diversity (Barreiro et al. 2009)
    'XPEHH_12_snp',  # XPEHH, max over two comparisons, (Sabeti et al. 2007, Grossman et al. 2010)
    'XPEHH_23_snp',  # XPEHH, max over two comparisons, (Sabeti et al. 2007, Grossman et al. 2010)
    'XPEHH_31_snp',  # XPEHH, max over two comparisons, (Sabeti et al. 2007, Grossman et al. 2010)
]

def check_epsilon(value, epsilon=0.0001, output=0.0001):
    """If less than epsilon, return epsilon"""
    if value < epsilon:
        return output
    else:
        return value
check_epsilon = np.vectorize(check_epsilon)

# not used
def get_FST_WC(n_1, n_2, p_1, p_2):
    """Return Weir and Cockerham's (1984) estimator of FST
        for bi-allelic SNPs, following Bhatia et al. (2013)"""
    a = 2*n_1*n_2/(n_1+n_2)
    b = (n_1*p_1*(1-p_1)+n_2*p_2*(1-p_2))/(n_1+n_2-2)
    c = (n_1*n_2/(n_1+n_2))*(p_1-p_2)**2
    FST_num = c-b
    FST_den = c-b+a*b
    # if negative or epsilon
    # FST_num = check_epsilon(FST_num, 0., 0.)  # round to 0 if neg
    # FST_den = check_epsilon(FST_den, 0., 0.)  # round to 0 if neg    
    # FST_snp = FST_num/check_epsilon(FST_den, epsilon, epsilon)
    # return FST components
    # return FST_snp, FST_num, FST_den
    return FST_num, FST_den

# not used
def get_FST_Nei(n_1, n_2, p_1, p_2):
    """Return Nei's ((1986) estimator of FST for bi-allelic 
        SNPs, following Bhatia et al. (2013)"""
    p_avg = (p_1+p_2)/2
    FST_num = (p_1-p_2)**2
    FST_den = (2*p_avg*(1-p_avg))
    # if negative or epsilon
    # FST_num = check_epsilon(FST_num, 0., 0.)  # round to 0 if neg
    # FST_den = check_epsilon(FST_den, 0., 0.)  # round to 0 if neg    
    # FST_snp = FST_num/check_epsilon(FST_den, epsilon, epsilon)
    # return FST components
    # return FST_snp, FST_num, FST_den
    return FST_num, FST_den

def get_FST_Hud(n_1, n_2, p_1, p_2):
    """Return Hudson et al.'s (1992) estimator of FST for 
        bi-allelic SNPs, following Bhatia et al. (2013)"""
    FST_num = (p_1-p_2)**2-p_1*(1-p_1)/(n_1-1)-p_2*(1-p_2)/(n_2-1)
    FST_den = p_1*(1-p_2)+p_2*(1-p_1)
    # if negative or epsilon
    # FST_num = check_epsilon(FST_num, 0., 0.)  # round to 0 if neg
    # FST_den = check_epsilon(FST_den, 0., 0.)  # round to 0 if neg    
    # FST_snp = FST_num/check_epsilon(FST_den, epsilon, epsilon)
    # return FST components
    # return FST_snp, FST_num, FST_den
    return FST_num, FST_den

alpha_list = dict()  # memoization
def get_alpha(n):
    """Fu (1995)"""
    # save small list of precomputed values in memory
    if n not in alpha_list:
        alpha_list[n] = np.sum([1/i for i in range(1, n)])
    return alpha_list[n]

beta_list = dict()  # memoization
def get_beta(n, i):
    """Fu (1995)"""
    assert 0 < i and i < n
    # save small list of precomputed values in memory
    if (n, i) not in beta_list:
        beta_list[(n, i)] = 2*n*(get_alpha(n+1)-get_alpha(i))/((n-i+1)*(n-i))-2/(n-i)
    return beta_list[(n, i)]

def get_sigma(n, i, j):
    """Fu (1995)"""
    assert 0 < i and i < n and 0 < j and j < n
    if i == j:
        if i < n/2:
            return get_beta(n,i+1)
        elif i == n/2:
            return 2*(get_alpha(n)-get_alpha(i))/(n-i)-1/(i*i)
        elif i > n/2:
            return get_beta(n,i)-1/(i*i)
        else:
            assert False
    else:
        if i < j:
            i, j = j, i
        if i+j < n:
            return (get_beta(n,i+1)-get_beta(n,i))/2
        elif i+j == n:
            return (get_alpha(n)-get_alpha(i))/(n-i)+(get_alpha(n)-
                get_alpha(j))/(n-j)-(get_beta(n,i)+get_beta(n,j+1))/2-1/(i*j) 
        elif i+j > n:
            return (get_beta(n,j)-get_beta(n,j+1))/2-1/(i*j)
        else:
            assert False

def get_const(theta_1_wts, theta_2_wts, n):
    """Return parameters for normalization given weights,
        following Achaz (2009) and Fu (1995)"""
    omega = theta_1_wts - theta_2_wts
    alpha = np.sum(np.arange(1,n)*np.square(omega))
    beta = np.sum([[i*j*omega[i-1]*omega[j-1]*get_sigma(n, i, j) \
        for i in np.arange(1, n)] for j in np.arange(1, n)])
    return alpha, beta

def get_SFS_norm(theta_1, theta_2, theta_S, alpha, beta):
    """Return neutrality test normalization given weights,
        following Achaz (2009) and Fu (1995)"""
    SFS_num = theta_1-theta_2
    SFS_den = np.sqrt(theta_S*(alpha+beta*theta_S))
    if SFS_den <= 0:
        return 'NA'
    else:
        return SFS_num/SFS_den

DIND_list = dict()  # memoization
def get_SFS_DIND(hap_X, filter_low, sfs_stats=True):
    DAC_window = np.sum(hap_X, axis=0)
    n, _ = hap_X.shape
    # get unfolded site frequency spectrum
    SFS_window = np.bincount(DAC_window)
    SFS_window = np.hstack((SFS_window, np.zeros(n+1-len(SFS_window))))[1:-1]
    SFS_window, seg_win = SFS_window*np.arange(1,n), sum(SFS_window)
    # theta weights
    if n not in DIND_list:
        theta_S_wts = np.array([1/i for i in range(1,n)])
        theta_Pi_wts = np.array([n-i for i in range(1,n)])
        theta_L_wts = np.array([1 for i in range(1,n)])
        theta_H_wts = np.array([i for i in range(1,n)])
        # # change weights for filter_low_freq_var
        # for i in range(filter_low):
        #     theta_S_wts[i] = 0
        #     theta_Pi_wts[i] = 0
        #     theta_L_wts[i] = 0
        #     theta_H_wts[i] = 0
        # normalize weights
        theta_S_wts = theta_S_wts/sum(theta_S_wts)
        theta_Pi_wts = theta_Pi_wts/sum(theta_Pi_wts)
        theta_L_wts = theta_L_wts/sum(theta_L_wts)
        theta_H_wts = theta_H_wts/sum(theta_H_wts)
        # SFS normalizations 
        if sfs_stats:
            a_Taj_D, b_Taj_D = get_const(theta_Pi_wts, theta_S_wts, n)
            a_FayWu_H, b_FayWu_H = get_const(theta_Pi_wts, theta_H_wts, n)
            a_Zeng_E, b_Zeng_E = get_const(theta_L_wts, theta_S_wts, n)
        # not needed for DIND
        else:
            a_Taj_D, a_FayWu_H, a_Zeng_E = my_arg_list('NA', 3)
            b_Taj_D, b_FayWu_H, b_Zeng_E = my_arg_list('NA', 3)
        # store in DIND table
        DIND_list[n] = [theta_S_wts, theta_Pi_wts, theta_L_wts, theta_H_wts, \
            a_Taj_D, a_FayWu_H, a_Zeng_E, b_Taj_D, b_FayWu_H, b_Zeng_E]
    else:
        [theta_S_wts, theta_Pi_wts, theta_L_wts, theta_H_wts, \
            a_Taj_D, a_FayWu_H, a_Zeng_E, b_Taj_D, b_FayWu_H, b_Zeng_E] = DIND_list[n]
    # get theta estimators
    theta_S = np.dot(SFS_window, theta_S_wts)
    theta_Pi = np.dot(SFS_window, theta_Pi_wts)
    theta_L = np.dot(SFS_window, theta_L_wts)
    theta_H = np.dot(SFS_window, theta_H_wts)
    return [theta_S, theta_Pi, theta_L, theta_H, seg_win, 
        theta_S_wts, theta_Pi_wts, theta_L_wts, theta_H_wts, 
        a_Taj_D, a_FayWu_H, a_Zeng_E, b_Taj_D, b_FayWu_H, b_Zeng_E]

def get_SFS_Delta(self, p, i, s):
    """Computes delta for thetas, mean log of ratios, 
        i.e. -log(theta_D/sqrt(theta_A1*theta_A2))"""
    s_1 = self.features_values[s][p][i]
    s_2 = self.features_values[s][(p+1)%3][i]
    s_3 = self.features_values[s][(p+2)%3][i]
    stat_list = [s_1, s_2, s_3]
    # case missing value
    if 'NA' in stat_list:
        return 'NA'
    else:
        # check epsilon before taking log
        # return np.log(s_3/check_epsilon(s_1)), np.log(s_2/check_epsilon(s_1))
        if s_1>0 and s_2>0:
            d_2 = np.log(s_2)-np.log(s_1)
        else:
            d_2 = 'NA'
        if s_1>0 and s_3>0:
            d_3 = np.log(s_3)-np.log(s_1)
        else:
            d_3 = 'NA'
        return d_2, d_3
        # return s_3/check_epsilon(s_1), s_2/check_epsilon(s_1)
        # # T_1, T_2, T_3 = -np.log(check_epsilon([s_1, s_2, s_3]))
        # # return T_1-T_2, T_1-T_3, T_1-(T_2+T_3)/2

def get_LD_pair(p_ij, p_i, p_j):
    """Return r^2_ij linkage disequilibrium"""
    return (p_ij-p_i*p_j)**2/(p_i*(1-p_i)*p_j*(1-p_j))

def get_rec_parallel(self, down_chunk, seed):
    np.random.seed(seed)
    features_values = self.features_values
    downsample_list = self.downsample_list
    pop_i = self.pop_i
    hap = self.data_shared['haplotypes'][self.pop_i]
    DAF = np.sum(hap, axis=0)/hap.shape[0]
    for i, d in zip(down_chunk, downsample_list[down_chunk]):
        low_i, high_i = self.get_window(d)
        recrates_win = np.mean(self.data_shared['recrates'][low_i:high_i])
        segsites_win = np.sum([np.abs(DAF[low_i:high_i]-0.5) < 0.5][0])
        features_values['recrates_win'][pop_i][i] = recrates_win
        features_values['segsites_win'][pop_i][i] = segsites_win
    return self

def get_freq_serial(self):
    features_values = self.features_values
    downsample_list = self.downsample_list
    pop_i = self.pop_i
    pop_list = self.pop_list
    # get DAF
    DAF_dict = dict()
    for p in self.pop_index:
        hap = self.data_shared['haplotypes'][p]
        n = hap.shape[0]
        DAF_dict[p] = np.sum(hap[:, :], axis=0)/n
    # # get DDAF_win
    # DDAF_list = DAF_dict[pop_i]-(DAF_dict[(pop_i-1)%3]+DAF_dict[(pop_i-2)%3])/2
    # for i, d in enumerate(downsample_list):
    #     low_i, high_i = self.get_window(d)
    #     # features_values['DDAF_win'][pop_i][i] = np.nanmax(DDAF_list[low_i:high_i])
    #     DDAF_win = np.count_nonzero(np.abs(DDAF_list[low_i:high_i]) > 0.3)/np.count_nonzero(
    #         ~np.isnan(DDAF_list[low_i:high_i]))  # # # 
    #     features_values['DDAF_win'][pop_i][i] = DDAF_win
    # get snp-level FSTs
    FST_num, FST_den = dict(), dict()
    samples = self.samples_list
    for p in self.pop_index:
        n = len(DAF_dict[p])
        FST_num[p], FST_den[p] = my_arg_list(my_arg_list('NA', n), 2)
        # get pairwise FST_snp, key: 0 with 01, 1 with 12, 2 with 20
        for i in range(len(DAF_dict[p])):
            FST_num[p][i], FST_den[p][i] = get_FST_Hud(
                samples[p], samples[(p+1)%3], DAF_dict[p][i], DAF_dict[(p+1)%3][i])
    self.FST_num, self.FST_den = FST_num, FST_den
    # FST_snp, FST_num, FST_den = dict(), dict(), dict()
    # samples = self.samples_list
    # for p in self.pop_index:
    #     n = len(DAF_dict[p])
    #     FST_snp[p], FST_num[p], FST_den[p] = my_arg_list(my_arg_list('NA', n), 3)
    #     # get pairwise FST_snp, key: 0 with 01, 1 with 12, 2 with 20
    #     for i in range(len(DAF_dict[p])):
    #         FST_snp[p][i], FST_num[p][i], FST_den[p][i] = get_FST_Hud(
    #             samples[p], samples[(p+1)%3], DAF_dict[p][i], DAF_dict[(p+1)%3][i])
    # self.FST_snp, self.FST_num, self.FST_den = FST_snp, FST_num, FST_den
    self.DAF_dict = DAF_dict
    return self

def get_freq_serial_standardize(self):
    features_values = self.features_values
    downsample_list = self.downsample_list
    pop_i = self.pop_i
    # get DAF
    DAF_dict = dict()
    for p in self.pop_index:
        hap = self.data_shared['haplotypes'][p]
        n = hap.shape[0]
        DAF_dict[p] = np.sum(hap[:, :], axis=0)/n
        features_values['DAF_snp'][p] = DAF_dict[p][downsample_list]
    # # get DDAF_win
    # DDAF_list = DAF_dict[pop_i]-(DAF_dict[(pop_i-1)%3]+DAF_dict[(pop_i-2)%3])/2
    # features_values['DDAF_snp'][pop_i] = DDAF_list[downsample_list]
    return self

def get_freq_parallel(self, down_chunk, seed):
    """Return allele frequency, fixation index, and genetic distance statistics:
        DAF, DDAF_win, FST_snp, FST_win, LSBL_win, PBS_win"""
    np.random.seed(seed)
    features_values = self.features_values
    downsample_list = self.downsample_list
    pop_i = self.pop_i
    # FST_snp, FST_num, FST_den = self.FST_snp, self.FST_num, self.FST_den
    FST_num, FST_den = self.FST_num, self.FST_den
    # get win-level FSTs
    # ascertaiment scheme: use SNPs polymorphic in any of the three populations
    for i, d in zip(down_chunk, downsample_list[down_chunk]):
        low_i, high_i = self.get_window(d)
        FST_win = dict()
        # get pairwise FST_win, key: 0 with 01, 1 with 12, 2 with 20
        for p in self.pop_index:
            # ratio of averages
            FST_win_num = np.sum(np.array(FST_num[p][low_i:high_i]))
            FST_win_den = np.sum(np.array(FST_den[p][low_i:high_i]))            
            # FST_win[p] = FST_win_num/check_epsilon(FST_win_den)
            if FST_win_den == 0:
                FST_win[p] = 'NA'
            else:
                FST_win[p] = FST_win_num/FST_win_den
                if FST_win[p] <= 0:
                    FST_win[p] = 0
                # if FST_win[p] >= 1:
                #     FST_win[p] = 1
        # get FST_win, LSBL_win, and PBS_win
        FST_win_12, FST_win_31, FST_win_23 = FST_win[pop_i], FST_win[(pop_i-1)%3], FST_win[(pop_i+1)%3]
        # get FST statistic
        # features_values['FST_win'][pop_i][i] = (FST_win_12 + FST_win_31)/2
        if FST_win_12 != 'NA' and FST_win_31 != 'NA':
            # print FST_win_12, FST_win_31
            features_values['FST_win'][pop_i][i] = (FST_win_12 + FST_win_31)/2
        # else:
        #     features_values['FST_win'][pop_i][i] = 'NA'
        # # get PBS statistic
        # T_12, T_31, T_23 = -np.log(check_epsilon(1-np.array([FST_win_12, FST_win_31, FST_win_23])))
        # PBS_win = check_epsilon((T_12+T_31-T_23)/2, 0., 0.)  # note, negative F3 a test for admixture
        # features_values['PBS_win'][pop_i][i] = PBS_win
    return self

def get_SFS_parallel(self, down_chunk, seed):
    """Return site frequency spectrum statistics:
        theta_Xi, theta_S, theta_Pi, theta_L, theta_H, 
        Taj_D, FuLi_F, FuLi_D, FayWu_H, Zeng_E, R2"""
    features_values = self.features_values
    downsample_list = self.downsample_list
    filter_low = self.params_shared.filter_low_freq_var
    pop_i = self.pop_i
    # get sfs features
    for p in self.pop_index:
        hap = self.data_shared['haplotypes'][p]
        n = hap.shape[0]
        # restrict to haplotype number >= 2 (else no SFS)
        if n >= 2:
            for i, d in zip(down_chunk, downsample_list[down_chunk]):
                low_i, high_i = self.get_window(d)
                hap_X = hap[:, low_i:high_i]
                DAC_window = np.sum(hap_X, axis=0)
                [theta_S, theta_Pi, theta_L, theta_H, seg_win, theta_S_wts, theta_Pi_wts, theta_L_wts, theta_H_wts, 
                    a_Taj_D, a_FayWu_H, a_Zeng_E, b_Taj_D, b_FayWu_H, b_Zeng_E] = get_SFS_DIND(hap_X, filter_low, True)
                # features_values['theta_S_win'][p][i] = theta_S
                features_values['theta_Pi_win'][p][i] = theta_Pi
                # # get SFS statistics
                if p == pop_i:
                    # if zero segregating sites then NA, if zero alpha and beta coeffs then NA;
                    # use windows only with at least two segregating sites, instead of one
                    # get SFS statistics
                    if seg_win >= 2:
                        Taj_D = get_SFS_norm(theta_Pi, theta_S, theta_S, a_Taj_D, b_Taj_D)
                        FayWu_H = get_SFS_norm(theta_Pi, theta_H, theta_S, a_FayWu_H, b_FayWu_H)
                        Zeng_E = get_SFS_norm(theta_L, theta_S, theta_S, a_Zeng_E, b_Zeng_E)
                        features_values['Taj_D_win'][p][i] = Taj_D
                        features_values['FayWu_H_win'][p][i] = FayWu_H
                        features_values['Zeng_E_win'][p][i] = Zeng_E
    # # get DIND-variant statistics
    # hap = self.data_shared['haplotypes'][self.pop_i]
    # n = hap.shape[0]
    # set_snp = set()
    # DIND_snp = dict()
    # for i, d in zip(down_chunk, downsample_list[down_chunk]):
    #     low_i, high_i = self.get_window(d)
    #     set_snp.update(range(low_i, high_i))
    # set_snp = np.sort(list(set_snp))
    # for i, d in zip(set_snp, self.data_shared['physpos'][set_snp]):
    #     low_i, high_i = self.get_window(d)
    #     hap_A = hap[np.nonzero(1-hap[:, i])[0], low_i:high_i]
    #     hap_D = hap[np.nonzero(hap[:, i])[0], low_i:high_i]
    #     A, D = hap_A.shape[0], hap_D.shape[0]
    #     # restrict to MAF >= 0.05 (same as iHS)
    #     MAF_min = 0.05
    #     DIND_snp[i] = np.nan
    #     if A/n >= MAF_min and D/n >= MAF_min:
    #         # restrict to haplotype number >= 2 (else no SFS)
    #         # if A >= 2 and D >= 2:
    #         if A >= 4 and D >= 4:
    #             # get derived and ancestral thetas
    #             theta_Pi_A = get_SFS_DIND(hap_A, filter_low, False)[1]
    #             theta_Pi_D = get_SFS_DIND(hap_D, filter_low, False)[1]
    #             if theta_Pi_D != 0:
    #                 # # # 
    #                 DIND_snp[i] = theta_Pi_A/theta_Pi_D
    #                 # # # 
    #             # else:
    #                 # DIND_snp[i] = 20
    #                 # # # 
    # for i, d in zip(down_chunk, downsample_list[down_chunk]):
    #     low_i, high_i = self.get_window(d)
    #     DIND_temp = [DIND_snp[j] for j in range(low_i, high_i)]
    #     if not np.isnan(DIND_temp).all():
    #         DIND_Pi_min = np.nanmin(DIND_temp)
    #         # # # 
    #         features_values['DIND_Pi_min'][pop_i][i] = DIND_Pi_min
    # get delta theta estimators
    for i, d in zip(down_chunk, downsample_list[down_chunk]):
        # D_theta_S_win_12, D_theta_S_win_31 = get_SFS_Delta(self, pop_i, i, 'theta_S_win')
        # features_values['D_theta_S_win'][pop_i][i] = max(D_theta_S_win_12, D_theta_S_win_31)
        D_theta_Pi_win_12, D_theta_Pi_win_31 = get_SFS_Delta(self, pop_i, i, 'theta_Pi_win')
        if D_theta_Pi_win_12 != 'NA' and D_theta_Pi_win_31 != 'NA':
            # print D_theta_Pi_win_12, D_theta_Pi_win_31
            features_values['D_theta_Pi_win'][pop_i][i] = (D_theta_Pi_win_12 + D_theta_Pi_win_31)/2
        # else:
        #     features_values['D_theta_Pi_win'][pop_i][i] = 'NA'
    return self

def get_hap_parallel(self, down_chunk, seed):
    """Return haplotype distribution statistics:
        K, avgHaplo, H1, H2, H12, H2/H1"""
    np.random.seed(seed)
    downsample_list = self.downsample_list
    features_values = self.features_values
    pop_i = self.pop_i
    hap = self.data_shared['haplotypes'][self.pop_i]
    # get hap features
    # # # get 1-HAF spectrum
    # # DAF = np.sum(hap, axis=0)
    for i, d in zip(down_chunk, downsample_list[down_chunk]):
        low_i, high_i = self.get_window(d)
        hap_temp = hap[:, low_i:high_i]
        # get haplotype freq distribution
        hap_dist = np.sort(np.unique([np.sum(np.left_shift(
            h, np.arange(len(h)))) for h in hap_temp], 
            return_counts=True)[1])[::-1]
        hap_dist = hap_dist/hap.shape[0]
        # get squared haplotype frequencies
        hap_square = np.square(hap_dist)
        # get K, H1, H2, H12, H2/H1
        K = float(len(hap_dist))
        H1 = sum(hap_square)
        if K >= 2:
            H2 = H1-hap_square[0]
            H12 = H1+(2*hap_dist[0]*hap_dist[1])
            H2H1 = H2/H1
        else:
            H2 = 0.
            H12 = H1
            H2H1 = 0.
        features_values['K_win'][pop_i][i] = K
        # # features_values['H1_win'][pop_i][i] = H1
        features_values['H2_win'][pop_i][i] = H2
        features_values['H12_win'][pop_i][i] = H12
        features_values['H2H1_win'][pop_i][i] = H2H1
        # # # get 1-HAF spectrum
        # # DAF_temp = DAF[low_i:high_i]
        # # HAF_1 = np.sort(np.dot(hap_temp, DAF_temp))[::-1]
    return self

def get_LD_parallel(self, down_chunk, seed):
    """Return linkage disequilibrium statistics:
        Wall_B, Wall_Q, ZA"""
    np.random.seed(seed)
    features_values = self.features_values
    downsample_list = self.downsample_list
    pop_i = self.pop_i
    hap = self.data_shared['haplotypes'][self.pop_i]
    # get ld features
    DAF = np.sum(hap, axis=0)/hap.shape[0]
    # get ZA, Wall_B, Wall_Q
    LD_list = dict()  # memoization
    for k, d in zip(down_chunk, downsample_list[down_chunk]):
        # keep only segsites
        low_i, high_i = self.get_window(d)
        window_seg = np.arange(low_i, high_i)[np.abs(DAF[low_i:high_i]-0.5) < 0.5]
        seg_win = len(window_seg)
        # need at least two segsites
        if seg_win >= 2:
            # get ZA
            ZA = 0
            for index, i in enumerate(window_seg[:-1]):
                j = window_seg[index+1]
                if i not in LD_list:
                    p_i = DAF[i]
                    p_j = DAF[j]
                    p_ij = np.sum(hap[:, i]*hap[:, j])/hap.shape[0]
                    LD_list[i] = get_LD_pair(p_ij, p_i, p_j)
                ZA += LD_list[i]
            ZA = ZA/(seg_win-1)
            features_values['ZA_win'][pop_i][k] = ZA
            # get Wall_B and Wall_Q
            B_prime = 0
            # A_partition = set()
            # adjacent segsites
            for ind, i in enumerate(window_seg[:-1]):
                j = window_seg[ind+1]
                dummy = 2*hap[:, i]+hap[:, j]
                # check congruency
                if len(np.unique(dummy)) == 2:
                    B_prime += 1
                    # binarize partition
                    dummy = 1*(dummy != dummy[0])
                    # # representation of partition
                    # # A_partition.add(np.sum(np.left_shift(dummy, np.arange(len(dummy)))))
            # save Wall_B and Wall_Q
            Wall_B = B_prime/(seg_win-1)
            features_values['Wall_B_win'][pop_i][k] = Wall_B
            # # Wall_Q = (B_prime+len(A_partition))/seg_win
            # # features_values['Wall_Q_win'][pop_i][k] = Wall_Q
    return self

def get_hapbin_serial(self, hapbin=False):
    """Return selscan integrated haplotype scores:
        iHS, D_iHH, nSL, D_nSL,
        XPEHH_12, XPEHH_23, XPEHH_31, XPEHH_final,
        D_XPEHH_12, D_XPEHH_23, D_XPEHH_31, D_XPEHH_final"""
    features_values = self.features_values
    downsample_list = self.downsample_list
    pop_i = self.pop_i
    # create posfile and hapfiles
    #   for simulations
    if not self.params_shared.isapplication:
        # selscan
        if hapbin:
            in_folder = '../../data/simulations/hapbin_in_temp/'
            out_folder = '../../data/simulations/hapbin_out_temp/'
        else:
            in_folder = '../../data/simulations/selscan_in_temp/'
            out_folder = '../../data/simulations/selscan_out_temp/'
        in_temp = '{0}_{1}'.format(self.meta_shared['label'], self.meta_shared['job'])
        self.meta_shared['pos_file'] = in_folder + in_temp + '.pos'
        self.meta_shared['hap_files'] = dict()
        for p_i in self.pop_index:
            self.meta_shared['hap_files'][p_i] = in_folder + in_temp + '.hap' + str(p_i+1)
        with open(self.meta_shared['pos_file'], 'w') as pos_file:
            chrm = 0  # placeholder for chrm
            for i, physpos in enumerate(self.data_shared['physpos']):
                varid = self.data_shared['varid'][i]
                genpos = self.data_shared['genpos'][i]
                pos_file.write(' '.join([str(chrm), str(varid), str(genpos), str(physpos)])+'\n')
        for p_i in self.pop_index:
            if hapbin:
                hap = self.data_shared['haplotypes'][p_i].transpose()
            else:
                hap = self.data_shared['haplotypes'][p_i]
            with open(self.meta_shared['hap_files'][p_i], 'w') as hap_file:
                hap_file.write('\n'.join([' '.join([str(i) for i in line]) for line in hap]))
    #   for applications, hapbin version
    else:
        if hapbin:
            in_folder = '../../data/applications/1000GP_Phase1_post/'
            out_folder = '../../data/applications/hapbin_out_temp/'
        else:
            return self
    # temporary file variables
    pos_file = self.meta_shared['pos_file']
    hap_file = self.meta_shared['hap_files'][pop_i]
    if not self.params_shared.isapplication:
        out_file = out_folder + '{0}_{1}'.format(self.meta_shared['label'], self.meta_shared['job'])
    else:
        out_file = out_folder + '{0}_{1}_{2}_{3}'.format(
            self.meta_shared['label'], self.meta_shared['job'], self.meta_shared['pop'], self.meta_shared['chr'])
    max_extend = int(self.params_shared.selscan_window_size*1e6/2)
    selscan_threads = self.params_shared.selscan_threads
    # compute iHS, D_iHH
    if hapbin:
        hapbin_iHS = 'time ihsbin --hap {0} --map {1} --out {2} '.format(hap_file, pos_file, out_file+'.ihs.out') + \
            '--minmaf 0.05 --cutoff 0.05 --max-extend {0} --binom'.format(max_extend)
        command_iHS = hapbin_iHS
    else:
        selscan_iHS = 'time selscan --ihs --hap {0} --map {1} --out {2} --maf 0.05 '.format(hap_file, pos_file, out_file) + \
            '--cutoff 0.05 --max-extend {0} --threads {1} --trunc-ok'.format(max_extend, selscan_threads)  # --keep-low-freq 
        command_iHS = selscan_iHS
    print ''
    print command_iHS
    os.system(command_iHS)
    # get iHS, D_iHH
    set_snp = set()
    for i, d in enumerate(downsample_list):
        low_i, high_i = self.get_window(d)
        set_snp.update(range(low_i, high_i))
    set_snp = np.sort(list(set_snp))
    var_list = dict(zip(self.data_shared['varid'][set_snp], set_snp))
    phys_list = dict(zip(self.data_shared['physpos'][set_snp], set_snp))
    iHS_snp_dict = dict(zip(set_snp, my_arg_list(np.nan, len(set_snp))))
    # D_iHH_snp_dict = dict(zip(set_snp, my_arg_list(np.nan, len(set_snp))))
    if not self.params_shared.isapplication:
        my_table_iHH = pd.read_csv("../simulations/standardize/my_simulations_iHH_{0}.txt".format(self.pop_list[pop_i]), sep='\t')
    else:
        my_table_iHH = pd.read_csv("../simulations/standardize/my_applications_iHH_{0}.txt".format(self.pop_list[pop_i]), sep='\t')
    with open(out_file + '.ihs.out', 'r') as f:
        if hapbin:
            for line in f.readlines()[1:]:
                line_list = line.split()
                varid = str(line_list[0])
                if varid in var_list:
                    DAF_snp = self.DAF_dict[pop_i][var_list[varid]]
                    iHH1, iHH0 = float(line_list[2]), float(line_list[1])  # iHH0 followed by iHH1
                    iHS_snp = np.log(iHH0)-np.log(iHH1)
                    # iHS_snp, D_iHH_snp = np.log(iHH0)-np.log(iHH1), abs(iHH0-iHH1)
                    iHS_snp = standardize_iHS_snp(iHS_snp, DAF_snp, self.pop_list[pop_i], my_table_iHH)
                    # D_iHH_snp = standardize_D_iHH_snp(D_iHH_snp, DAF_snp, self.pop_list[pop_i], my_table_iHH)
                    iHS_snp_dict[var_list[varid]] = iHS_snp
                    # D_iHH_snp_dict[var_list[varid]] = D_iHH_snp
        else:
            for line in f.readlines():
                line_list = line.split()
                physpos = int(line_list[1])
                if physpos in phys_list:
                    DAF_snp = self.DAF_dict[pop_i][phys_list[physpos]]
                    iHH1, iHH0 = float(line_list[3]), float(line_list[4])
                    iHS_snp = np.log(iHH0)-np.log(iHH1)
                    # iHS_snp, D_iHH_snp = np.log(iHH0)-np.log(iHH1), abs(iHH0-iHH1)
                    iHS_snp = standardize_iHS_snp(iHS_snp, DAF_snp, self.pop_list[pop_i], my_table_iHH)
                    # D_iHH_snp = standardize_D_iHH_snp(D_iHH_snp, DAF_snp, self.pop_list[pop_i], my_table_iHH)
                    iHS_snp_dict[phys_list[physpos]] = iHS_snp
                    # D_iHH_snp_dict[phys_list[physpos]] = D_iHH_snp
    for i, d in enumerate(downsample_list):
        low_i, high_i = self.get_window(d)
        iHS_list = np.array([iHS_snp_dict[j] for j in range(low_i, high_i)])
        # D_iHH_list = np.array([D_iHH_snp_dict[j] for j in range(low_i, high_i)])
        iHS_list = iHS_list[np.isfinite(iHS_list)]
        # D_iHH_list = D_iHH_list[np.isfinite(D_iHH_list)]
        if len(iHS_list) > 0:
            iHS_win = np.count_nonzero(np.abs(iHS_list) > 2)/len(iHS_list)  # # # 
            features_values['iHS_win'][pop_i][i] = iHS_win
        # if len(D_iHH_list) > 0:
        #     D_iHH_win = np.count_nonzero(np.abs(D_iHH_list) > 2)/len(D_iHH_list)  # # # 
        #     features_values['D_iHH_win'][pop_i][i] = D_iHH_win

    # compute XPEHH, D_XPEHH
    pop_i_temp = [(self.pop_i, (self.pop_i+1)%3), ((self.pop_i+2)%3, self.pop_i)]
    for q1, q2 in pop_i_temp:
        hap_file = self.meta_shared['hap_files'][q1]
        ref_file = self.meta_shared['hap_files'][q2]
        if not self.params_shared.isapplication:
            out_file = out_folder + '{0}_{1}_{2}{3}'.format(self.meta_shared['label'], self.meta_shared['job'], q1+1, q2+1)
        else:
            out_file = out_folder + '{0}_{1}_{2}_{3}{4}'.format(self.meta_shared['label'], self.meta_shared['pop'], self.meta_shared['chr'], q1+1, q2+1)
        if hapbin:
            hapbin_XPEHH = 'time xpehhbin --hapA {0} --hapB {1} --map {2} --out {3} '.format(hap_file, ref_file, pos_file, out_file+'.xpehh.out') + \
                '--cutoff 0.05 --max-extend {0} --binom'.format(max_extend, self.params_shared.selscan_threads)  # --wagh  # # #  
            command_XPEHH = hapbin_XPEHH
        else:
            selscan_XPEHH = 'time selscan --xpehh --hap {0} --ref {1} --map {2} --out {3} '.format(hap_file, ref_file, pos_file, out_file) + \
                '--cutoff 0.05 --max-extend {0} --threads {1} --trunc-ok'.format(max_extend, self.params_shared.selscan_threads)  # --wagh  --keep-low-freq  # # #  
            command_XPEHH = selscan_XPEHH
        print ''
        print command_XPEHH
        os.system(command_XPEHH)
    # get XPEHH, D_XPEHH
    XPEHH_dict = dict()
    if not self.params_shared.isapplication:
        my_table_XPEHH = pd.read_csv("../simulations/standardize/my_simulations_XPEHH_{0}.txt".format(self.pop_list[pop_i]), sep='\t')
    else:
        my_table_XPEHH = pd.read_csv("../simulations/standardize/my_applications_XPEHH_{0}.txt".format(self.pop_list[pop_i]), sep='\t')
    for q1, q2 in pop_i_temp:
        XPEHH_dict[(q1, q2)] = dict(zip(set_snp, my_arg_list(np.nan, len(set_snp))))
        if not self.params_shared.isapplication:
            out_file = out_folder + '{0}_{1}_{2}{3}'.format(self.meta_shared['label'], self.meta_shared['job'], q1+1, q2+1)
        else:
            out_file = out_folder + '{0}_{1}_{2}_{3}{4}'.format(self.meta_shared['label'], self.meta_shared['pop'], self.meta_shared['chr'], q1+1, q2+1)
        with open(out_file + '.xpehh.out', 'r') as f:
            if hapbin:
                for line in f.readlines()[1:]:
                    line_list = line.split()
                    varid = str(line_list[0])
                    if varid in var_list:
                        iHHA, iHHB = float(line_list[1]), float(line_list[2])  # iHHA followed by iHHB
                        XPEHH_snp = np.log(iHHA)-np.log(iHHB) 
                        XPEHH_snp = standardize_XPEHH_snp(XPEHH_snp, self.pop_list[pop_i], q1, q2, my_table_XPEHH)
                        XPEHH_dict[(q1, q2)][var_list[varid]] = XPEHH_snp
            else:
                for line in f.readlines()[1:]:
                    line_list = line.split()
                    physpos = int(line_list[1])
                    if physpos in phys_list:
                        iHHA, iHHB = float(line_list[4]), float(line_list[6])
                        XPEHH_snp = np.log(iHHA)-np.log(iHHB)
                        XPEHH_snp = standardize_XPEHH_snp(XPEHH_snp, self.pop_list[pop_i], q1, q2, my_table_XPEHH)
                        XPEHH_dict[(q1, q2)][phys_list[physpos]] = XPEHH_snp
    # find windowed-average across sites
    for i, d in enumerate(downsample_list):
        low_i, high_i = self.get_window(d)
        XPEHH_12_list = np.array([XPEHH_dict[(self.pop_i, (self.pop_i+1)%3)][j] for j in range(low_i, high_i)])
        XPEHH_31_list = np.array([XPEHH_dict[((self.pop_i+2)%3, self.pop_i)][j] for j in range(low_i, high_i)])
        XPEHH_12_list = XPEHH_12_list[np.isfinite(XPEHH_12_list)]
        XPEHH_31_list = XPEHH_31_list[np.isfinite(XPEHH_31_list)]
        if len(XPEHH_12_list > 0) and len(XPEHH_31_list) > 0:
            XPEHH_12_win = np.count_nonzero(np.abs(XPEHH_12_list) > 2.5)/len(XPEHH_12_list)
            XPEHH_31_win = np.count_nonzero(np.abs(XPEHH_31_list) > 2.5)/len(XPEHH_31_list)
            features_values['XPEHH_win'][pop_i][i] = (XPEHH_12_win + XPEHH_31_win)/2

            # 
            # 
            # 

    return self

def get_hapbin_serial_standardize(self, hapbin=False):
    """Return selscan integrated haplotype scores:
        iHS, D_iHH, nSL, D_nSL,
        XPEHH_12, XPEHH_23, XPEHH_31, XPEHH_final,
        D_XPEHH_12, D_XPEHH_23, D_XPEHH_31, D_XPEHH_final"""
    features_values = self.features_values
    downsample_list = self.downsample_list
    pop_i = self.pop_i
    # create posfile and hapfiles
    #   for simulations
    if not self.params_shared.isapplication:
        # selscan
        if hapbin:
            in_folder = '../../data/simulations/hapbin_in_temp/'
            out_folder = '../../data/simulations/hapbin_out_temp/'
        else:
            in_folder = '../../data/simulations/selscan_in_temp/'
            out_folder = '../../data/simulations/selscan_out_temp/'
        in_temp = '{0}_{1}'.format(self.meta_shared['label'], self.meta_shared['job'])
        self.meta_shared['pos_file'] = in_folder + in_temp + '.pos'
        self.meta_shared['hap_files'] = dict()
        for p_i in self.pop_index:
            self.meta_shared['hap_files'][p_i] = in_folder + in_temp + '.hap' + str(p_i+1)
        with open(self.meta_shared['pos_file'], 'w') as pos_file:
            chrm = 0  # placeholder for chrm
            for i, physpos in enumerate(self.data_shared['physpos']):
                varid = self.data_shared['varid'][i]
                genpos = self.data_shared['genpos'][i]
                pos_file.write(' '.join([str(chrm), str(varid), str(genpos), str(physpos)])+'\n')
        for p_i in self.pop_index:
            if hapbin:
                hap = self.data_shared['haplotypes'][p_i].transpose()
            else:
                hap = self.data_shared['haplotypes'][p_i]
            with open(self.meta_shared['hap_files'][p_i], 'w') as hap_file:
                hap_file.write('\n'.join([' '.join([str(i) for i in line]) for line in hap]))
    #   for applications, hapbin version
    else:
        if hapbin:
            in_folder = '../../data/applications/1000GP_Phase1_post/'
            out_folder = '../../data/applications/hapbin_out_temp/'
        else:
            return self
    # temporary file variables
    pos_file = self.meta_shared['pos_file']
    hap_file = self.meta_shared['hap_files'][pop_i]
    if not self.params_shared.isapplication:
        out_file = out_folder + '{0}_{1}'.format(self.meta_shared['label'], self.meta_shared['job'])
    else:
        out_file = out_folder + '{0}_{1}_{2}_{3}'.format(
            self.meta_shared['label'], self.meta_shared['job'], self.meta_shared['pop'], self.meta_shared['chr'])
    max_extend = int(self.params_shared.selscan_window_size*1e6/2)
    selscan_threads = self.params_shared.selscan_threads
    # compute iHS, D_iHH
    if hapbin:
        hapbin_iHS = 'time ihsbin --hap {0} --map {1} --out {2} '.format(hap_file, pos_file, out_file+'.ihs.out') + \
            '--minmaf 0.05 --cutoff 0.05 --max-extend {0} --binom'.format(max_extend)
        command_iHS = hapbin_iHS
    else:
        selscan_iHS = 'time selscan --ihs --hap {0} --map {1} --out {2} --maf 0.05 '.format(hap_file, pos_file, out_file) + \
            '--cutoff 0.05 --max-extend {0} --threads {1} --trunc-ok'.format(max_extend, selscan_threads)  # --keep-low-freq 
        command_iHS = selscan_iHS
    print ''
    print command_iHS
    os.system(command_iHS)
    # get iHS, D_iHH
    var_list = dict(zip(self.data_shared['varid'][downsample_list], range(len(downsample_list))))
    phys_list = dict(zip(self.data_shared['physpos'][downsample_list], range(len(downsample_list))))
    with open(out_file + '.ihs.out', 'r') as f:
        if hapbin:
            for line in f.readlines()[1:]:
                line_list = line.split()
                varid = str(line_list[0])
                if varid in var_list:
                    iHH1, iHH0 = float(line_list[2]), float(line_list[1])  # iHH0 followed by iHH1
                    iHS_snp = np.log(iHH0)-np.log(iHH1)
                    # iHS_snp, D_iHH_snp = np.log(iHH0)-np.log(iHH1), abs(iHH0-iHH1)
                    features_values['iHS_snp'][pop_i][var_list[varid]] = iHS_snp
                    # features_values['D_iHH_snp'][pop_i][var_list[varid]] = D_iHH_snp
        else:
            for line in f.readlines():
                line_list = line.split()
                physpos = int(line_list[1])
                if physpos in phys_list:
                    # DAF_snp = float(line_list[2])
                    iHH1, iHH0 = float(line_list[3]), float(line_list[4])
                    iHS_snp = np.log(iHH0)-np.log(iHH1)
                    # iHS_snp, D_iHH_snp = np.log(iHH0)-np.log(iHH1), abs(iHH0-iHH1)
                    features_values['iHS_snp'][pop_i][phys_list[physpos]] = iHS_snp
                    # features_values['D_iHH_snp'][pop_i][phys_list[physpos]] = D_iHH_snp
    # compute XPEHH, D_XPEHH
    pop_i_temp = [(self.pop_i, (self.pop_i+1)%3), ((self.pop_i+2)%3, self.pop_i)]
    for q1, q2 in pop_i_temp:
        hap_file = self.meta_shared['hap_files'][q1]
        ref_file = self.meta_shared['hap_files'][q2]
        if not self.params_shared.isapplication:
            out_file = out_folder + '{0}_{1}_{2}{3}'.format(self.meta_shared['label'], self.meta_shared['job'], q1+1, q2+1)
        else:
            out_file = out_folder + '{0}_{1}_{2}_{3}{4}'.format(self.meta_shared['label'], self.meta_shared['pop'], self.meta_shared['chr'], q1+1, q2+1)
        if hapbin:
            hapbin_XPEHH = 'time xpehhbin --hapA {0} --hapB {1} --map {2} --out {3} '.format(hap_file, ref_file, pos_file, out_file+'.xpehh.out') + \
                '--cutoff 0.05 --max-extend {0} --binom'.format(max_extend, self.params_shared.selscan_threads)  # --wagh  # # #  
            command_XPEHH = hapbin_XPEHH
        else:
            selscan_XPEHH = 'time selscan --xpehh --hap {0} --ref {1} --map {2} --out {3} '.format(hap_file, ref_file, pos_file, out_file) + \
                '--cutoff 0.05 --max-extend {0} --threads {1} --trunc-ok'.format(max_extend, self.params_shared.selscan_threads)  # --wagh  --keep-low-freq  # # #  
            command_XPEHH = selscan_XPEHH
        print ''
        print command_XPEHH
        os.system(command_XPEHH)
    # get XPEHH, D_XPEHH
    for q1, q2 in pop_i_temp:
        if not self.params_shared.isapplication:
            out_file = out_folder + '{0}_{1}_{2}{3}'.format(self.meta_shared['label'], self.meta_shared['job'], q1+1, q2+1)
        else:
            out_file = out_folder + '{0}_{1}_{2}_{3}{4}'.format(self.meta_shared['label'], self.meta_shared['pop'], self.meta_shared['chr'], q1+1, q2+1)
        with open(out_file + '.xpehh.out', 'r') as f:
            if hapbin:
                for line in f.readlines()[1:]:
                    line_list = line.split()
                    varid = str(line_list[0])
                    if varid in var_list:
                        iHHA, iHHB = float(line_list[1]), float(line_list[2])  # iHHA followed by iHHB
                        XPEHH_snp = np.log(iHHA)-np.log(iHHB)
                        features_values['XPEHH_{0}{1}_snp'.format(q1+1, q2+1)][pop_i][var_list[varid]] = XPEHH_snp
            else:
                for line in f.readlines()[1:]:
                    line_list = line.split()
                    physpos = int(line_list[1])
                    if physpos in phys_list:
                        iHHA, iHHB = float(line_list[4]), float(line_list[6])
                        XPEHH_snp = np.log(iHHA)-np.log(iHHB)
                        features_values['XPEHH_{0}{1}_snp'.format(q1+1, q2+1)][pop_i][phys_list[physpos]] = XPEHH_snp
    return self

def get_parallel(self, get_worker, num_proc):
    """Partitions downsample and calls given worker on each partition"""
    downsample_subs = np.array_split(range(len(self.downsample_list)), num_proc)
    p_list = dict()
    for i, down_chunk in enumerate(downsample_subs):
        seed = np.random.random_integers(0, np.iinfo(np.int32).max)
        p_list[i] = multiprocessing.Process(target=get_worker, args=(self, down_chunk, seed))
        p_list[i].start()
    for i, down_chunk in enumerate(downsample_subs):
        p_list[i].join()
    return self

class Features(object):
    """Class for calculating, storing, and printing summary statistics"""
    def __init__(self, params_shared, data_shared, meta_shared):
        self.meta_shared = meta_shared
        self.data_shared = data_shared
        self.params_shared = params_shared
        if not self.params_shared.isapplication:
            self.params_sims = self.params_shared
        self.samples_list = [
            np.shape(self.data_shared['haplotypes'][0])[0], 
            np.shape(self.data_shared['haplotypes'][1])[0], 
            np.shape(self.data_shared['haplotypes'][2])[0]]
        self.pop = self.meta_shared['pop']
        self.pop_list = self.params_shared.pop_list
        self.pop_i = self.pop_list.index(self.pop)
        self.pop_index = range(len(self.pop_list))

    def get_window(self, i):
        """Helper method to get indices for window statistics"""
        # get window indices
        low = i - 1e6/2*self.params_shared.window_stat_size
        high = i + 1e6/2*self.params_shared.window_stat_size
        low_i = np.searchsorted(self.data_shared['physpos'], low, side='left')
        high_i = np.searchsorted(self.data_shared['physpos'], high, side='right')
        return (low_i, high_i)

    def get_downsample(self):
        """Downsample sites, ignoring chromosome margins"""
        # margin indices
        if self.params_shared.isapplication:
            low = self.data_shared['physpos'][0] + 1e6*self.params_shared.chromosome_margin
            high = self.data_shared['physpos'][-1] - 1e6*self.params_shared.chromosome_margin
        else:
            low = 1e6*self.params_shared.chromosome_margin
            high = 1e6*self.params_shared.chromosome_length - low
        # use win step-size
        if self.params_shared.isapplication:
            window_step_size = self.params_shared.window_step_size_app
        else:
            window_step_size = self.params_shared.window_step_size_sim
        self.downsample_list = np.arange(low+1e6*window_step_size/2, 
            high-1e6*window_step_size/2+0.1, 1e6*window_step_size)
        # check valid window
        downsample_temp = []
        for i, d in enumerate(self.downsample_list):
            low_i, high_i = self.get_window(d)
            if low_i < high_i:
                downsample_temp.append(d)
            else:
                print i, d, low_i, high_i
        self.downsample_list = np.array(downsample_temp)
        return self

    def get_downsample_standardize(self):
        """Downsample sites, ignoring chromosome margins"""
        # margin indices
        if self.params_shared.isapplication:
            low = self.data_shared['physpos'][0] + 1e6*self.params_shared.chromosome_margin
            high = self.data_shared['physpos'][-1] - 1e6*self.params_shared.chromosome_margin
        else:
            low = 1e6*self.params_shared.chromosome_margin
            high = 1e6*self.params_shared.chromosome_length - low
        # use win step-size
        low_i = np.argmax(self.data_shared['physpos'] >= low)
        high_i = np.argmax(self.data_shared['physpos'] > high)
        self.downsample_list = np.arange(low_i, high_i)
        return self

    def get_init_features(self):
        # initializing
        self.features_values = dict()  # shared
        features_values = self.features_values
        downsample_list = self.downsample_list
        pop_i = self.pop_i
        # for metadata
        for s in features_list:
            features_values[s] = dict()
            features_values[s][pop_i] = my_arg_list('NA', len(downsample_list))
        # for features
        for s in features_stat:
            if s in ['theta_Xi_win', 'theta_S_win', 'theta_Pi_win', 'theta_L_win', 'theta_H_win']:
                features_values[s] = my_arg_list(None, len(self.pop_index))
                for p in self.pop_index:
                    features_values[s][p] = multiprocessing.Array('d', my_arg_list(my_NA, len(downsample_list)))
            else:
                features_values[s] = dict()
                features_values[s][pop_i] = multiprocessing.Array('d', my_arg_list(my_NA, len(downsample_list)))
        # for p in self.pop_index:
        for i, d in enumerate(downsample_list):
            physpos_d = d
            features_values['key'][pop_i][i] = self.meta_shared['key']
            features_values['physpos'][pop_i][i] = physpos_d
            # # features_values['varid'][pop_i][i] = self.data_shared['varid'][d]
            # # features_values['physpos'][pop_i][i] = self.data_shared['physpos'][d]
            # # features_values['genpos'][pop_i][i] = self.data_shared['genpos'][d]
            # # features_values['recsite'][pop_i][i] = self.data_shared['recrates'][d]
            if not self.params_shared.isapplication:
                features_values['dist_physpos'][pop_i][i] = self.params_sims.physpos_s - physpos_d
                # # physpos_d = self.data_shared['physpos'][d]
                # # genpos_d = self.data_shared['genpos'][d]
                # # features_values['dist_genpos'][pop_i][i] = abs(self.params_sims.genpos_s - genpos_d)
                # # features_values['dist_physpos'][pop_i][i] = abs(self.params_sims.physpos_s - physpos_d)
                if self.params_sims.stype == 'neutral':
                    site_class = 'neutral'
                else:
                    if self.params_sims.stype == 'hardsweep':
                        site_class = 'hardlinked'
                        if self.params_sims.physpos_s == physpos_d:
                            site_class = 'hardsweep'
                    if self.params_sims.stype == 'softsweep':
                        site_class = 'softlinked'
                        if self.params_sims.physpos_s == physpos_d:
                            site_class = 'softsweep'
                features_values['site_class'][pop_i][i] = site_class
        return self

    def get_init_features_standardize(self):
        # initializing
        self.features_values = dict()  # shared
        features_values = self.features_values
        downsample_list = self.downsample_list
        pop_i = self.pop_i
        # for metadata
        for s in features_list:
            features_values[s] = dict()
            features_values[s][pop_i] = my_arg_list('NA', len(downsample_list))
        # for features
        for s in features_standardize:
            if s in ['theta_Xi_win', 'theta_S_win', 'theta_Pi_win', 'theta_L_win', 'theta_H_win']:
                features_values[s] = my_arg_list(None, len(self.pop_index))
                for p in self.pop_index:
                    features_values[s][p] = multiprocessing.Array('d', my_arg_list(my_NA, len(downsample_list)))
            else:
                features_values[s] = dict()
                features_values[s][pop_i] = multiprocessing.Array('d', my_arg_list(my_NA, len(downsample_list)))
        # for p in self.pop_index:
        for i, d in enumerate(downsample_list):
            # physpos_d = d
            features_values['key'][pop_i][i] = self.meta_shared['key']
            # # features_values['physpos'][pop_i][i] = physpos_d
            # # features_values['varid'][pop_i][i] = self.data_shared['varid'][d]
            features_values['physpos'][pop_i][i] = self.data_shared['physpos'][d]
            # # features_values['genpos'][pop_i][i] = self.data_shared['genpos'][d]
            # # features_values['recsite'][pop_i][i] = self.data_shared['recrates'][d]
            if not self.params_shared.isapplication:
                physpos_d = self.data_shared['physpos'][d]
                features_values['dist_physpos'][pop_i][i] = abs(self.params_sims.physpos_s - physpos_d)
                # # physpos_d = self.data_shared['physpos'][d]
                # # genpos_d = self.data_shared['genpos'][d]
                # # features_values['dist_genpos'][pop_i][i] = abs(self.params_sims.genpos_s - genpos_d)
                # # features_values['dist_physpos'][pop_i][i] = abs(self.params_sims.physpos_s - physpos_d)
                if self.params_sims.stype == 'neutral':
                    site_class = 'neutral'
                else:
                    if self.params_sims.stype == 'hardsweep':
                        site_class = 'hardlinked'
                        if self.params_sims.physpos_s == physpos_d:
                            site_class = 'hardsweep'
                    if self.params_sims.stype == 'softsweep':
                        site_class = 'softlinked'
                        if self.params_sims.physpos_s == physpos_d:
                            site_class = 'softsweep'
                features_values['site_class'][pop_i][i] = site_class
        return self

    def get_final_features(self, hapbin):
        print ''
        num_proc = 1  # 4  # ~X-speedup, ~X-memory  # # # 
        print 'rec-p'
        get_parallel(self, get_rec_parallel, num_proc)
        print 'freq-s'
        get_freq_serial(self)  # must do first
        print 'freq-p'
        get_parallel(self, get_freq_parallel, num_proc)
        print 'sfs-p'
        get_parallel(self, get_SFS_parallel, num_proc)
        print 'hap-p'
        get_parallel(self, get_hap_parallel, num_proc)
        print 'ld-p'
        get_parallel(self, get_LD_parallel, num_proc)
        print 'ehh-s'
        get_hapbin_serial(self, hapbin)
        return self

    def get_final_features_standardize(self, hapbin):
        print ''
        num_proc = 1  # 4  # ~X-speedup, ~X-memory  # # # 
        print 'freq-s'
        get_freq_serial_standardize(self)  # must do first
        print 'ehh-s'
        get_hapbin_serial_standardize(self, hapbin)
        return self

class Metadata(object):
    """Class for calculating, storing, and printing simulation metadata"""
    def __init__(self, params_shared, meta_shared):
        self.meta_shared = meta_shared
        self.params_shared = params_shared
        if not self.params_shared.isapplication:
            self.params_sims = self.params_shared

    def get_final_metadata(self):
        self.metadata_values = dict()
        metadata_values = self.metadata_values
        for s in metadata_list:
            metadata_values[s] = 'NA'
        if self.params_shared.isapplication:
            metadata_values['key'] = self.meta_shared['key']
            metadata_values['pop'] = self.meta_shared['pop']
            metadata_values['chr'] = self.meta_shared['chr']
        else:
            metadata_values['key'] = self.params_sims.key
            metadata_values['pop'] = self.params_sims.pop
            metadata_values['stype'] = self.params_sims.stype
            metadata_values['model'] = self.params_sims.model
            metadata_values['seed'] = self.params_sims.seed
            metadata_values['mutrate'] = self.params_sims.mutrate
            metadata_values['recrate'] = self.params_sims.recrate
            metadata_values['spop'] = self.params_sims.spop
            metadata_values['scoeff'] = self.params_sims.scoeff
            metadata_values['sfreq'] = self.params_sims.sfreq
            metadata_values['stime'] = time_y(self.params_sims.stime, self.params_sims.Nref, self.params_sims.Tgen)
            metadata_values['complete'] = self.params_sims.complete
            metadata_values['sweep_freq_bin'] = self.params_sims.sweep_freq_bin_list[self.params_sims.pop_i]
            metadata_values['sweep_freq_pop'] = self.params_sims.sweep_freq_pop_list[self.params_sims.pop_i]
            metadata_values['sweep_fix_time'] = self.params_sims.sweep_fix_time_list[self.params_sims.pop_i]
            metadata_values['sweep_dur_time'] = self.params_sims.sweep_dur_time_list[self.params_sims.pop_i]
            metadata_values['origins'] = self.params_sims.origins
        return self

def get_features(params_shared, data_shared, meta_shared, hapbin=True):
    """Output features for given empirical data"""
    return Features(params_shared, data_shared, meta_shared).get_downsample().get_init_features().get_final_features(hapbin=hapbin)

def get_standardize(params_shared, data_shared, meta_shared, hapbin=True):
    """Output training SNPs for standardization"""
    return Features(params_shared, data_shared, meta_shared \
        ).get_downsample_standardize().get_init_features_standardize().get_final_features_standardize(hapbin=hapbin)

def get_metadata(params_shared, meta_shared):
    """Output features for given empirical data"""
    return Metadata(params_shared, meta_shared).get_final_metadata()
