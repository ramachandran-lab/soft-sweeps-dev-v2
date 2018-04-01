#------------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University           #
# Created 2016-01-08, Modified 2017-10-25                    #
# Conversions between DADI model units and MSMS model units  #
#------------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
from copy import copy
import numpy as np
import pandas as pd

def msms_seed():
    """Generate seed input for MSMS"""
    a, b = np.random.random_integers(np.iinfo(np.int32).min, np.iinfo(np.int32).max, size=2)
    return str(a+(b<<32))

def uniform_sample(minv, maxv):
    """Samples uniformly between min and max of prior."""
    return np.random.uniform(minv, maxv)

def uniform_logarithmic(minv, maxv):
    """Samples log-uniformly between min and max of prior"""
    if minv > 0 and maxv > 0:
        return np.exp(np.random.uniform(np.log(minv), np.log(maxv)))
    elif minv == 0 and maxv == 0:
        return 0
    else:
        assert False

def uniform_sample_discrete(minv, maxv, steps):
    """Samples uniformly between min and max of prior."""
    if steps == None:
        return uniform_sample(minv, maxv)
    step_size = (maxv-minv)/steps
    return np.random.choice(np.arange(minv, maxv+step_size, step_size))

def uniform_logarithmic_discrete(minv, maxv, steps):
    """Samples log-uniformly between min and max of prior"""
    if steps == None:
        return uniform_logarithmic(minv, maxv)
    step_size = (np.log(maxv)-np.log(minv))/steps
    if minv > 0 and maxv > 0:
        return np.exp(np.random.choice(np.arange(
            np.log(minv), np.log(maxv)+step_size, step_size)))
    elif minv == 0 and maxv == 0:
        return 0
    else:
        assert False

def bitruncated_exponential(mean, minv, maxv):
    """Sample from exponential distribution between min and max values"""
    temp = -1
    while temp < minv or temp > maxv:
        temp = np.random.exponential(mean)
    return temp

def round_str(number, decimals=0):
    """Round float to int if not NA"""
    if number != 'NA':
        if decimals != 0:
            return "{0:.{1}f}".format(number, decimals)
        else:
            return int(number)
    return 'NA'

def time_4N(time, Nref, Tgen):
    """Convert time from years to units of 4N gen"""
    if time != 'NA':
        return time/(4*Nref*Tgen)
    return 'NA'
time_4N = np.vectorize(time_4N)

def time_y(time, Nref, Tgen):
    """Convert time from units of 4N gen to years"""
    if time != 'NA':
        return 4*Nref*Tgen*time
    return 'NA'
time_y = np.vectorize(time_y)

def coeff_2Ns(coeff, Nref):
    """Convert coeff from units of s to 2Ns"""
    if coeff != 'NA':
        return 2*Nref*coeff
    return 'NA'
coeff_2Ns = np.vectorize(coeff_2Ns)

def coeff_s(coeff, Nref):
    """Convert coeff from units of 2Ns to s"""
    if coeff != 'NA':
        return coeff/(2*Nref)
    return 'NA'
coeff_s = np.vectorize(coeff_s)

def size_Nref(Neff, Nref):
    """Convert Neff from ratio of Nref to ind"""
    if Neff != 'NA':
        return Neff/(Nref)
    return 'NA'
size_Nref = np.vectorize(size_Nref)

def size_ind(Neff, Nref):
    """Convert Neff from ind to ratio of Nref"""
    if Neff != 'NA':
        return Neff*Nref
    return 'NA'
size_ind = np.vectorize(size_ind)

def migr_4Nm(migr, Nref):
    """Convert migr from units of m to 4Nm"""
    if migr != 'NA':
        return 4*Nref*migr
    return 'NA'
migr_4Nm = np.vectorize(migr_4Nm)

def migr_m(migr, Nref):
    """Convert migr from units of 4Nm to m"""
    if migr != 'NA':
        return migr/(4*Nref)
    return 'NA'
migr_m = np.vectorize(migr_m)

def rate_4Nx(rate, chrmlen, Nref):
    """Convert from rate per bp to units of theta or rho = 4Nx,
        where x is mu for mutrate, x is r for recrate"""
    if rate != 'NA':
        return 1e6*4*Nref*rate*chrmlen
    return 'NA'
rate_4Nx = np.vectorize(rate_4Nx)

def rate_bp(rate, chrmlen, Nref):
    """Convert from units of theta or rho to rate per bp = 4Nx,
        where x is mu for mutrate, x is r for recrate"""
    if rate != 'NA':
        return rate/(1e6*4*Nref*chrmlen)
    return 'NA'
rate_bp = np.vectorize(rate_bp)

def rate_cMMb(rate, chrmlen, Nref):
    """Convert from units of rho to cM/Mb,
        where x is r for recrate"""
    if rate != 'NA':
        return rate/(4*Nref*0.01*chrmlen)
    return 'NA'

def growth_exp(growth, Nref):
    """Convert percent per generation to exp growth rate in 4N gen"""
    if growth != 'NA':
        return 4*Nref*np.log(1+0.01*growth)
    return 'NA'
growth_exp = np.vectorize(growth_exp)

def growth_per(growth, Nref):
    """Convert exp growth rate in 4N gen to percent per generation"""
    if growth != 'NA':
        return 100*(np.exp(growth/(4*Nref))-1)
    return 'NA'
growth_per = np.vectorize(growth_per)

def size_final(Ninit, growth, time):
    """Calculates Neff final from Neff initial as ratio of Nref,
        exp growth rate, and elapsed (forward) time in 4N gen units"""
    if Ninit != 'NA' and growth != 'NA':
        return Ninit*np.exp(growth*time)
    return 'NA'
size_final = np.vectorize(size_final)

def size_init(Nfinal, growth, time):
    """Calculates Neff initial from Neff final as ratio of Nref,
        exp growth rate, and elapsed (backward) time in 4N gen units"""
    if Nfinal != 'NA' and growth != 'NA':
        return Nfinal*np.exp(-growth*time)
    return 'NA'
size_init = np.vectorize(size_init)

def get_physical_pos(msms_position, chrmlen):
    """Converts msms float position to int physical position"""
    return int(np.ceil(1e6*msms_position*chrmlen))
get_physical_pos = np.vectorize(get_physical_pos)

def get_genetic_pos_old(msms_position, recrate, Nref):
    """Converts msms float position to genetic position,
        under assumption of constant recombination rate"""
    return 1e2*msms_position*recrate/(4*Nref)
get_genetic_pos_old = np.vectorize(get_genetic_pos_old)

def get_genetic_pos_new(physical_pos, recrate, chrmlen, Nref):
    """Converts physical position to genetic position,
        under assumption of constant recombination rate per simulation,
        recall that recombination rate varies across simulations"""
    return 1e-4*physical_pos*recrate/(4*chrmlen*Nref)
get_genetic_pos_new = np.vectorize(get_genetic_pos_new)

def round_to_n(x, n=1):
    """Rounds x to n sig figs"""
    if not np.isclose(x, 0.):
        return np.round(x, -np.int(np.floor(np.log10(abs(x))-n+1)))
    else:
        return 0.
round_to_n = np.vectorize(round_to_n)

def my_arg_list(value, n):
    """Returns a list of value of length n"""
    return [copy(value) for i in range(n)]

my_NA = -1000000000.
def check_main(value, sfigs, stat):
    """Use instead of 'NA' for multiprocessing arrays for floats"""
    if isinstance(value, float):
        if np.isclose(value, my_NA):
            return 'NA'
        # Do not round position information
        elif stat in ["physpos", "dist_physpos"]:
            return str(value)
        else:
            return str(round_to_n(value, sfigs))
    else:
        return str(value)

def standardize_iHS_snp(iHS_value, DAF_value, pop, my_table):
    DAF_bins = int(np.ceil(40*DAF_value))
    if np.isclose(DAF_value, 0.05):  # include 0.05 in 2nd bin
        DAF_bins += 1
    iHS_snp_mean = my_table.loc[my_table['DAF_bins']==DAF_bins][['iHS_snp_mean']].as_matrix()[0][0]
    iHS_snp_sd = my_table.loc[my_table['DAF_bins']==DAF_bins][['iHS_snp_sd']].as_matrix()[0][0]
    return (iHS_value-iHS_snp_mean)/iHS_snp_sd

def standardize_D_iHH_snp(D_iHH_value, DAF_value, pop, my_table):
    DAF_bins = int(np.ceil(40*DAF_value))
    if np.isclose(DAF_value, 0.05):  # include 0.05 in 2nd bin
        DAF_bins += 1
    D_iHH_snp_mean = my_table.loc[my_table['DAF_bins']==DAF_bins][['D_iHH_snp_mean']].as_matrix()[0][0]
    D_iHH_snp_sd = my_table.loc[my_table['DAF_bins']==DAF_bins][['D_iHH_snp_sd']].as_matrix()[0][0]
    return (D_iHH_value-D_iHH_snp_mean)/D_iHH_snp_sd

def standardize_XPEHH_snp(XPEHH_value, pop, q1, q2, my_table):
    XPEHH_snp_mean = my_table[['XPEHH_{0}{1}_snp_mean'.format(q1+1, q2+1)]].as_matrix()[0][0]
    XPEHH_snp_sd = my_table[['XPEHH_{0}{1}_snp_sd'.format(q1+1, q2+1)]].as_matrix()[0][0]
    return (XPEHH_value-XPEHH_snp_mean)/XPEHH_snp_sd
