#---------------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University              #
# Created 2015-07-23, Modified 2017-11-13                       #
# Demography and selection parameters for training simulations  #
#---------------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import numpy as np

class ParametersDemography(object):
    def __init__(self):
        # model = 'GravelFinal2011'
        model = 'Jouganos2017-YCC'

        #  Gravel et al.'s (2011) exome model
        if model == 'GravelExome':
            self.model = 'GravelExome'  # description for demographic model being used
            self.Tgen = 25.  # years per generation
            self.Nref = 10000.  # reference population size
            self.NA = 0.731  # effective population sizes
            self.NAF = 1.5388
            self.NB = 0.2758
            self.NEU_f = 3.31934
            self.NAS_f = 2.62955
            self.rEU = 107.854  # population growth rates
            self.rAS = 123.808
            self.mAFB = 8.  # forward migration rates
            self.mAFEU = 0.68
            self.mAFAS = 0.232
            self.mEUAS = 2.36
            self.mBAF = 8.  # backward migration rates
            self.mEUAF = 0.68
            self.mASAF = 0.232
            self.mASEU = 2.36
            self.tAF = 0.316  # population event times
            self.tB = 0.098
            self.tEUAS = 0.028
            self.mutrate = 944.

        #  Gravel et al.'s (2011) intergenic model
        if model == 'GravelIntergenic':
            self.model = 'GravelIntergenic'  # description for demographic model being used
            self.Tgen = 25.  # years per generation
            self.Nref = 10000.  # reference population size
            self.NA = 0.73  # effective population sizes
            self.NAF = 1.23
            self.NB = 0.21
            self.NEU_f = 2.95249
            self.NAS_f = 5.34035
            self.rEU = 159.681  # population growth rates
            self.rAS = 219.397
            self.mAFB = 10.  # forward migration rates
            self.mAFEU = 1.2
            self.mAFAS = 0.76
            self.mEUAS = 3.84
            self.mBAF = 10.   # backward migration rates
            self.mEUAF = 1.2
            self.mASAF = 0.76
            self.mASEU = 3.84
            self.tAF = 0.22  # population event times
            self.tB = 0.14
            self.tEUAS = 0.0212
            self.mutrate = 944.

        #  Gravel et al.'s (2011) low-coverage + exome model
        if model == 'GravelFinal2011':
            self.model = 'GravelFinal2011'  # description for demographic model being used
            self.Tgen = 25.  # years per generation
            self.Nref = 10000.  # reference population size
            self.NA = 0.731  # effective population sizes
            self.NAF = 1.4474
            self.NB = 0.1861
            self.NEU_f = 3.38139
            self.NAS_f = 4.53697
            self.rEU = 151.712  # population growth rates
            self.rAS = 191.541
            self.mAFB = 6.  # forward migration rates
            self.mAFEU = 1.
            self.mAFAS = 0.312
            self.mEUAS = 1.244
            self.mBAF = 6.  # backward migration rates
            self.mEUAF = 1.
            self.mASAF = 0.312
            self.mASEU = 1.244
            self.tAF = 0.148  # population event times
            self.tB = 0.051
            self.tEUAS = 0.023
            self.mutrate = 944.

        #  Gravel et al.'s (2011) low-coverage + exome model
        if model == 'Jouganos2017-YCC':
            self.model = 'Jouganos2017-YCC'  # description for demographic model being used
            self.Tgen = 29.  # years per generation
            self.Nref = 10000.  # reference population size
            self.NA = 1.1273  # effective population sizes
            self.NAF = 2.3721
            self.NB = 0.3104
            self.NEU_f = 3.95007
            self.NAS_f = 8.31915
            self.rEU = 78.3233  # population growth rates
            self.rAS = 123.409
            self.mAFB = 6.32  # forward migration rates
            self.mAFEU = 0.44
            self.mAFAS = 0.192
            self.mEUAS = 1.676
            self.mBAF = 6.32  # backward migration rates
            self.mEUAF = 0.44
            self.mASAF = 0.192
            self.mASEU = 1.676
            self.tAF = 0.2690  # population event times
            self.tB = 0.1078
            self.tEUAS = 0.0365
            self.mutrate = 576.

        #  Gravel et al.'s (2011) low-coverage + exome model
        if model == 'Jouganos2017-YGC':
            self.model = 'Jouganos2017-YGC'  # description for demographic model being used
            self.Tgen = 29.  # years per generation
            self.Nref = 10000.  # reference population size
            self.NA = 1.1273  # effective population sizes
            self.NAF = 2.4486
            self.NB = 0.3034
            self.NEU_f = 3.40428
            self.NAS_f = 9.01976
            self.rEU = 67.9423  # population growth rates
            self.rAS = 119.820
            self.mAFB = 6.24  # forward migration rates
            self.mAFEU = 0.4
            self.mAFAS = 0.192
            self.mEUAS = 1.596
            self.mBAF = 6.24  # backward migration rates
            self.mEUAF = 0.4
            self.mASAF = 0.192
            self.mASEU = 1.596
            self.tAF = 0.3009  # population event times
            self.tB = 0.1043
            self.tEUAS = 0.0379
            self.mutrate = 576.

        #  Gravel et al.'s (2011) low-coverage + exome model
        if model == 'Jouganos2017-YGK':
            self.model = 'Jouganos2017-YGK'  # description for demographic model being used
            self.Tgen = 29.  # years per generation
            self.Nref = 10000.  # reference population size
            self.NA = 1.1273  # effective population sizes
            self.NAF = 2.3908
            self.NB = 0.2986
            self.NEU_f = 3.08403
            self.NAS_f = 8.07514
            self.rEU = 67.9423  # population growth rates
            self.rAS = 119.820
            self.mAFB = 6.4  # forward migration rates
            self.mAFEU = 0.408
            self.mAFAS = 0.24
            self.mEUAS = 1.784
            self.mBAF = 6.4  # backward migration rates
            self.mEUAF = 0.408
            self.mASAF = 0.24
            self.mASEU = 1.784
            self.tAF = 0.2698  # population event times
            self.tB = 0.1034
            self.tEUAS = 0.0371
            self.mutrate = 576.

        #  Gravel et al.'s (2011) low-coverage + exome model
        if model == 'Jouganos2017-LCC':
            self.model = 'Jouganos2017-LCC'  # description for demographic model being used
            self.Tgen = 29.  # years per generation
            self.Nref = 10000.  # reference population size
            self.NA = 1.1273  # effective population sizes
            self.NAF = 2.9034
            self.NB = 0.2746
            self.NEU_f = 3.86593
            self.NAS_f = 8.08216
            self.rEU = 83.9119  # population growth rates
            self.rAS = 131.783
            self.mAFB = 7.08  # forward migration rates
            self.mAFEU = 0.6
            self.mAFAS = 0.252
            self.mEUAS = 1.772
            self.mBAF = 7.08  # backward migration rates
            self.mEUAF = 0.6
            self.mASAF = 0.252
            self.mASEU = 1.772
            self.tAF = 0.1879  # population event times
            self.tB = 0.0853
            self.tEUAS = 0.0345
            self.mutrate = 576.

class ParametersShared(ParametersDemography):
    def __init__(self):
        super(ParametersShared, self).__init__()
        self.isapplication = True  # check if application data or simulation data

        # Feature calculations
        self.pop_list = ['YRI', 'CEU', 'CHB']
        self.samples_AF = 176  # 20  # haplotype sample size from AF  # # # 
        self.samples_EU = 170  # 20  # haplotype sample size from EU  # # # 
        self.samples_AS = 194  # 20  # haplotype sample size from AS  # # # 
        self.window_stat_size = 0.05  # size of windows for summary statistics in Mb
        self.chromosome_margin = 0.125  # ignore ends of chromosome, compare selscan_window_size #
        self.selscan_window_size = 0.20  # twice selscan param max-extend, selscan default = 2.0
        self.selscan_save_inputs = False  # save selscan files for each execution separately
        self.selscan_threads = 32  # selscan param for number of threads, selscan default = 1
        self.window_step_size_app = 0.025  # 0.01  # size of step-size for overlappping windows in Mb  # # # 
        self.filter_low_freq_var = 2  # 0 filter monomorphic, 1 filter singletons, 2 filter doubletons  # 

class ParametersSimulations(ParametersShared):
    def __init__(self):
        super(ParametersSimulations, self).__init__()
        self.isapplication = False
        self.remove_spos = True  # False  # # # 
        
        # Sampling parameters
        self.msms_threads = 32  # msms param for number of threads, msms default = 1
        self.simulations_ne = 50  # neutral simulations per population  # 
        self.simulations_hs = 50  # hard sweep simulations per population  # 
        self.simulations_ss = 50  # soft sweep simulations per population  # 
        self.window_step_size_stan = 0.05  # size of step-size for overlappping windows in Mb  # # # 
        self.window_step_size_train = 0.05  # size of step-size for overlappping windows in Mb  # # # 
        self.window_step_size_link = 0.025  # 0.01  # size of step-size for overlappping windows in Mb  # # # 
        self.chromosome_length_stan = 1.25  # size of our simulated genomic region in Mb  # # # 
        self.chromosome_length_train = 0.30  # size of our simulated genomic region in Mb  # # # 
        self.chromosome_length_link = 1.30  # size of our simulated genomic region in Mb  # # # 
        
        # Genomic parameters
        mutr = self.mutrate
        self.mutrate_prior = (mutr, mutr, None)  # uniform prior for theta per Mb  # 1.44e-8 (Gravel et al. 2013) # # # 
        self.recrate_prior = (440., 80., 8000.)  # log uniform prior for rho per Mb  # range 0.2 to 20 cM/Mb, mean 1.1e-8 # 
        
        # Selection parameters
        self.spos_prior = (0.5, 0.5, None)  # the msms location of sweep between 0 and 1, to 5 sig figs  #
        self.sdom_prior = (0.5, 0.5, None)  # constant prior for dominance coefficient h
        self.fixation_cutoff = 0.95  # frequency threshold for defining a completed sweep
        # self.low_sweep_cutoff = 0.05  # threshold on minimum final population frequency of sweep
        self.low_sweep_cutoff = 0.2  # threshold on minimum final population frequency of sweep
        self.sweepbin_discrete = 10  # round sweep population frequency to nearest tenth

        # Hard sweep parameters
        # self.hstime_prior = (0.001, 'tB', None) # earliest time for hard sweep, float, 'tB', or 'tEUAS'
        self.hstime_prior = (0.00431, 0.03017, None) # earliest time for hard sweep, float, 'tB', or 'tEUAS' # from 5 to 35kya
        self.hscoeff_prior = (0.0005, 0.5, None)  # [0.05, 0.5]  # log uniform prior for selection coefficient s  # 
        # self.hsfreq_prior = (0.00005, 0.05, None)  # log uniform prior for initial freq, (1/2*Nref)
        self.hsfreq_prior = (0.00005, 0.2, None)  # log uniform prior for initial freq, (1/2*Nref)

        # Soft sweep parameters
        # self.sstime_prior = (0.001, 'tB', None) # earliest time for hard sweep, float, 'tB', or 'tEUAS'
        self.sstime_prior = (0.00431, 0.03017, None) # earliest time for hard sweep, float, 'tB', or 'tEUAS' # from 5 to 35kya
        self.sscoeff_prior = (0.0005, 0.5, None)  # [0.05, 0.5]  # log uniform prior for selection coefficient s  # 
        # self.ssfreq_prior = (0.00005, 0.05, None)  # log uniform prior for initial freq, (1/2*Nref)
        self.hsfreq_prior = (0.00005, 0.2, None)  # log uniform prior for initial freq, (1/2*Nref)
