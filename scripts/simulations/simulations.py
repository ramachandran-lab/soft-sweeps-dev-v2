#--------------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University             #
# Created 2015-07-23, Modified 2017-10-25                      #
# Generates MSMS simulations from given parameters and priors  #
#--------------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import re, subprocess
import numpy as np
from parameters import *
from conversions import *

msms_prefix = 'time java -jar msms.jar '
sweep_suffix = ' -oOC -oTrace -SFC'  # -SFC conditions on allele not being lost

class NeutralModel(object):
    """Class for neutral models"""
    def __init__(self, params_sims=ParametersSimulations(), key='NA', pop='NA', seed='NA'):
        for attr in params_sims.__dict__:
            setattr(self, attr, getattr(params_sims, attr))
        self.key = key
        self.pop = pop
        self.seed = seed
        if seed == 'NA':
            self.seed = msms_seed()
        self.pop_i = self.pop_list.index(self.pop)
        # identifiers
        self.stype = 'neutral'
        # genomic rate parameters
        self.mutrate = round_to_n(uniform_sample_discrete(*self.mutrate_prior)*self.chromosome_length, 4)  # # # 
        self.recrate = round_to_n(bitruncated_exponential(*self.recrate_prior)*self.chromosome_length, 4)  # # # 
        # selection parameters
        self.spop = '111'
        self.stime = 'NA'
        self.scoeff = 'NA'
        self.sfreq = 'NA'
        self.sdom = 'NA'
        self.spos = 'NA'
        self.physpos_s = 'NA'
        self.genpos_s = 'NA'
        # attributes for MSMS command and output
        self.msms = 'NA'
        self.content = 'NA'
        # postprocess attributes for NeutralModel
        self.segsites = 'NA'
        self.msms_pos_list = 'NA'
        self.physpos_list = 'NA'
        self.genpos_list = 'NA'
        self.haplotypes_list = 'NA'
        # postprocess attributes for SweepModel
        self.origins = 0  # 'NA'
        self.complete = 'NA'
        self.frequency_list = 'NA'
        self.sweep_freq_pop_list = my_arg_list('NA', len(self.pop_list))
        self.sweep_freq_bin_list = my_arg_list('NA', len(self.pop_list))
        self.sweep_fix_time_list = my_arg_list('NA', len(self.pop_list))
        self.sweep_dur_time_list = my_arg_list('NA', len(self.pop_list))
        # errors and troubleshooting
        self.check = True
        self.problems = set()

    def construct_msms(self):
        """Construct MSMS command for NeutralModel object"""
        # overall parameters
        samples_total = self.samples_AF + self.samples_EU + self.samples_AS
        msms_1 = '-N {0} -ms {1} 1 -I 3 {2} {3} {4} '.format(
            self.Nref, samples_total, self.samples_AF, self.samples_EU, self.samples_AS)
        # rate parameters
        msms_2 = '-t {0} -r {1} '.format(self.mutrate, self.recrate)
        # initial growth and size parameters
        msms_3 = '-n 1 {0} -n 2 {1} -n 3 {2} -g 2 {3} -g 3 {4} '.format(
            self.NAF, self.NEU_f, self.NAS_f, self.rEU, self.rAS, )
        # initial migration parameters
        msms_4 = '-m 1 2 {0} -m 2 1 {1} -m 1 3 {2} -m 3 1 {3} -m 2 3 {4} -m 3 2 {5} '.format(
            self.mAFEU, self.mEUAF, self.mAFAS, self.mASAF, self.mEUAS, self.mASEU)
        # interm parameters
        msms_5 = '-ej {0} 3 2 -en {0} 2 {1} -em {0} 1 2 {2} -em {0} 2 1 {3} '.format(
            self.tEUAS, self.NB, self.mAFB, self.mBAF)
        # final parameters
        msms_6 = '-ej {0} 2 1 -en {1} 1 {2}'.format(self.tB, self.tAF, self.NA)
        # concatenate
        msms_neutral = msms_1 + msms_2 + msms_3 + msms_4 + msms_5 + msms_6
        # add MSMS seed
        if self.stype == 'neutral':
            msms_neutral += ' -seed ' + str(self.seed) + ' -threads ' + str(self.msms_threads)
        # add MSMS prefix
        self.msms = msms_prefix + msms_neutral
        return self

    def print_msms(self):
        """Print MSMS command to terminal"""
        # must construct before print
        self.construct_msms()
        print ''
        print self.msms
        return self

    def print_problems(self):
        print ''
        for i in sorted(self.problems):
            print i

    def execute_msms(self, outfile='NA'):
        """Execute MSMS command for NeutralModel object"""
        # must construct before execute
        self.construct_msms()
        # execute msms command and save output
        self.content = subprocess.Popen([self.msms], stdout=subprocess.PIPE, 
            shell=True).communicate()[0].split('\n')
        # save to file if location given
        if outfile != 'NA':
            with open(outfile, 'w') as f:
                for line in self.content:
                    f.write(line+'\n')
        return self

    def postprocess_msms(self, outfile='NA'):
        """Postprocess MSMS content for NeutralModel object"""
        # must execute before postprocess
        self.execute_msms(outfile)
        # postprocess msms content and save parts
        for i, line in enumerate(self.content):
            # segsites
            if re.search('segsites: ', line):
                self.segsites = int(re.split(': ', line)[1])
                # msms_pos_list
                self.msms_pos_list = np.array(
                    re.split(': ', self.content[i+1])[1].split(), dtype=np.float)
                # varid_list
                self.varid_list = np.array([str(s) for s in range(self.segsites)], dtype=object)
                # recrates_list
                self.recrates_list = np.array(my_arg_list(rate_cMMb(self.recrate, self.chromosome_length, self.Nref), 
                    self.segsites))
                # haplotypes_list
                self.haplotypes_list = []
                self.haplotypes_list.append(np.array([list(hapi) for hapi in self.content[
                    i+2:i+2+self.samples_AF]]))
                self.haplotypes_list.append(np.array([list(hapi) for hapi in self.content[
                    i+2+self.samples_AF:i+2+self.samples_AF+self.samples_EU]]))
                self.haplotypes_list.append(np.array([list(hapi) for hapi in self.content[
                    i+2+self.samples_AF+self.samples_EU:i+2+self.samples_AF+self.samples_EU+self.samples_AS]]))
                # pseudo spos for neutral
                if self.stype == 'neutral':
                    self.spos = round_to_n(uniform_sample_discrete(*self.spos_prior), 4)
                # find nearest position to spos
                self.spos_i = np.abs(self.msms_pos_list-self.spos).argmin()
                # check, due to rounding, spos may appear more than once in msms_pos_list
                spos_count = sum(self.msms_pos_list == self.msms_pos_list[self.spos_i])
                if spos_count != 1:
                    self.check = False
                    self.problems.add('- spos count of {0} is not equal to 1'.format(spos_count))
                # cleanup, get physical positions of snps,
                # note, since we rounded float representations to nearest integer base pair,
                # this can result in neighboring snps being assigned the same physical position,
                # we correct for this by shifting the snp on the right up by 1 physical position 
                # if after spos index, and left by 1 physical position if before spos index
                self.physpos_list = np.array(get_physical_pos(self.msms_pos_list, self.chromosome_length))
                phys_count = 0
                for j in np.where(self.msms_pos_list > self.spos)[0][:-1]:
                    if self.physpos_list[j] >= self.physpos_list[j+1]:
                        self.physpos_list[j+1] = self.physpos_list[j]+1
                        phys_count += 1
                for j in np.where(self.msms_pos_list < self.spos)[0][1:][::-1]:
                    if self.physpos_list[j] <= self.physpos_list[j-1]:
                        self.physpos_list[j-1] = self.physpos_list[j]-1
                        phys_count += 1                        
                self.problems.add('+ corrected {0} ({1}%) of physpos'.format(
                    phys_count, round_to_n(100*phys_count/self.segsites, 2)))
                # save physical location of spos and genetic location of spos
                self.genpos_list = np.array(get_genetic_pos_new(self.physpos_list, self.recrate, self.chromosome_length, 
                    self.Nref))
                self.physpos_s = 1e6*self.spos*self.chromosome_length
                self.genpos_s = get_genetic_pos_new(self.physpos_s, self.recrate, self.chromosome_length, self.Nref)
                # convert haplotypes from chars to ints
                if self.stype == 'neutral':
                    self.haplotypes_list = [happ.astype(int) for happ in self.haplotypes_list]
                # print problems with simulations
                if self.stype == 'neutral':
                    self.print_problems()
        return self

def construct_SI(spop, sfreq):
    """Returns -SI flag"""
    return ' '.join([str(sfreq*int(i)) for i in list(spop)])

def construct_Sc(spop, pop, scoeff, sdom, Nref):
    """Returns -Sc flag, after converting to MSMS units"""
    scoeff = coeff_2Ns(scoeff, Nref)
    if int(spop[pop-1]) == 1:
        return '{0} {1} 0'.format(2*scoeff, 2*scoeff*sdom)
    return '0 0 0'

class SweepModel(NeutralModel):
    """Class for selective sweep models"""
    def __init__(self, params_sims, stype, key='NA', pop='NA', seed='NA'):
        """Construct NeutralModel object with sampled selection parameters"""
        super(SweepModel, self).__init__(params_sims, key, pop, seed)
        # identifiers
        self.stype = stype
        # selection priors
        if self.stype == 'hardsweep':
            stime_lower, stime_upper, stime_discrete, scoeff_prior, sfreq_prior = \
            self.hstime_prior[0], self.hstime_prior[1], self.hstime_prior[2], self.hscoeff_prior, self.hsfreq_prior
        elif self.stype == 'softsweep':
            stime_lower, stime_upper, stime_discrete, scoeff_prior, sfreq_prior = \
            self.sstime_prior[0], self.sstime_prior[1], self.sstime_prior[2], self.sscoeff_prior, self.ssfreq_prior
        else:
            raise AssertionError, 'no stype assigned'
        # set earliest sweep
        if stime_upper == 'tB':
            stime_upper = self.tB
        elif stime_upper == 'tEUAS':
            stime_upper = self.tEUAS
        # do sanity checks
        if stime_upper <= stime_lower:
            raise AssertionError, 'stime_upper <= stime_lower'
        if stime_upper <= 0:
            raise AssertionError, 'stime_upper <= 0'
        if stime_lower <= 0:
            raise AssertionError, 'stime_lower <= 0'
        stime_prior = (stime_lower, stime_upper, stime_discrete)
        # stime, scoeff, sfreq, spos, sdom
        self.stime = round_to_n(uniform_sample_discrete(*stime_prior), 4)
        self.scoeff = round_to_n(uniform_logarithmic_discrete(*scoeff_prior), 4)
        self.sfreq = round_to_n(uniform_logarithmic_discrete(*sfreq_prior), 4)
        self.spos = round_to_n(uniform_sample_discrete(*self.spos_prior), 4)
        self.sdom = round_to_n(uniform_sample_discrete(*self.sdom_prior), 4)
        # spop parameter is then defined
        if self.pop == self.pop_list[0]:
            if self.stime < self.tB:
                self.spop = '100'
        elif self.pop == self.pop_list[1]:
            if self.stime < self.tEUAS:
                self.spop = '010'
            elif self.stime < self.tB:
                self.spop = '011'
        elif self.pop == self.pop_list[2]:
            if self.stime < self.tEUAS:
                self.spop = '001'
            elif self.stime < self.tB:
                self.spop = '011'
        
    def construct_msms(self):
        """Construct MSMS command for SweepModel object"""
        # neutral parameters
        msms_neutral = super(SweepModel, self).construct_msms().msms
        # sweep time, location, and initial frequency
        msms_sweep = ' -SI {0} 3 {1} -Sp {2} -Smark'.format(
            self.stime, construct_SI(self.spop, self.sfreq), self.spos)
        # sweep coefficients in each population
        msms_sweep += ' -Sc 0 1 {0} -Sc 0 2 {1} -Sc 0 3 {2}'.format(
            construct_Sc(self.spop, 1, self.scoeff, self.sdom, self.Nref), 
            construct_Sc(self.spop, 2, self.scoeff, self.sdom, self.Nref), 
            construct_Sc(self.spop, 3, self.scoeff, self.sdom, self.Nref))
        # concatenate, add seed and suffix
        self.msms = msms_neutral+msms_sweep+' -seed '+str(self.seed)+\
            ' -threads '+str(self.msms_threads)+sweep_suffix
        return self

    def postprocess_msms(self, outfile='NA'):
        """Postprocess MSMS content for SweepModel object"""
        # postprocess NeutralModel first
        super(SweepModel, self).postprocess_msms(outfile)
        # origins
        for line in self.content:
            if re.search('OriginCount:', line):
                self.origins = int(re.split(':', line)[1])
                break
        # frequency trajectory
        begin_i = None
        end_i = None
        for i, line in enumerate(self.content):
            if re.search('Frequency Trace:', line):
                begin_i = i+2
                break
        for i, line in enumerate(self.content):
            if re.search('segsites: ', line):
                end_i = i
                break
        self.frequency_list = np.array(
            [line.split() for line in self.content[begin_i:end_i]], dtype=float)
        # drop columns, for freq trajectory
        self.frequency_list = np.delete(self.frequency_list, [1, 3, 5], 1)
        # convert times, for freq trajectory
        self.frequency_list[:,0] = time_y(self.frequency_list[:,0], self.Nref, self.Tgen)
        # save population final frequency (not sample derived allele frequency)
        self.sweep_freq_pop_list = self.frequency_list[-1][1:]
        self.sweep_freq_bin_list = np.rint(self.sweepbin_discrete*self.sweep_freq_pop_list)/self.sweepbin_discrete
        # check if sweep is complete or incomplete (compare vs fixation_cutoff)
        if self.sweep_freq_pop_list[self.pop_i] > self.fixation_cutoff:
            self.complete = 1
        else:
            self.complete = 0
        # check, population derived allele frequency at sweep
        if all([i == 0 for i in self.sweep_freq_pop_list]):
            self.check = False
        if all([i == 1 for i in self.sweep_freq_pop_list]):
            self.check = False
        # pop freq of beneficial mutation
        n, m = self.frequency_list.shape
        # lower threshold for pop freq of sweep
        # must increase by low_sweep_cutoff compared to sfreq
        sweep_pop_freq = self.frequency_list[-1, 1+self.pop_i]
        # must increase by low_sweep_cutoff
        if self.low_sweep_cutoff is not None:
            if sweep_pop_freq < self.low_sweep_cutoff:
                self.check = False
                self.problems.add('- sweep {0} lower than {1}'.format(
                    round_to_n(sweep_pop_freq, 2), round_to_n(self.low_sweep_cutoff, 2)))
        # find and save fixation times of sweeps in past
        for j in range(1, m):
            if self.frequency_list[-1, j] < self.fixation_cutoff:
                continue
            # find earliest time freq > cutoff
            for i in range(n-1)[::-1]:
                if self.frequency_list[i, j] > self.fixation_cutoff:
                    self.sweep_fix_time_list[j-1] = self.frequency_list[i, 0]
            # fix pop 3 if before pop 2 and 3 split
            if self.sweep_fix_time_list[2] == time_y(self.tEUAS, self.Nref, self.Tgen):
                self.sweep_fix_time_list[2] = self.sweep_fix_time_list[1]
            # fix pop 2 if before pop 1 and 2 split
            if self.sweep_fix_time_list[1] == time_y(self.tB, self.Nref, self.Tgen):
                self.sweep_fix_time_list[1] = self.sweep_fix_time_list[0]
            # note: sweep_fix_time is earliest time when frequency_list > fixation_cutoff
            # but only if sweep_freq_pop (freq at present) is also > fixation_cutoff
            # this checks for sweeps that reach fixation_cutoff but drift below due to chance
            for i in range(len(self.spop)):
                if self.sweep_freq_pop_list[i] > self.fixation_cutoff:
                    continue
                else:  # revert back to 'NA'
                    self.sweep_fix_time_list[i] = 'NA'
            # find and save duration times of sweeps
            for i in range(len(self.spop)):
                if self.sweep_fix_time_list[i] != 'NA':
                    self.sweep_dur_time_list[i] = time_y(self.stime, self.Nref, self.Tgen)-self.sweep_fix_time_list[i]
        # check, that hard sweeps have 1 origin and soft sweeps have >1 origin
        if self.stype == 'hardsweep' and self.origins > 1:
            # self.check = False
            self.problems.add('+ hs origins of {0} is not equal to 1'.format(self.origins))
        elif self.stype == 'softsweep' and self.origins < 2:
            # self.check = False
            self.problems.add('+ ss origins of {0} is not less than 2'.format(self.origins))
        # cleanup, change non-unity characters in soft sweep haplotypes to unity
        # also, check that sweep site is not fixed in all pops and thus is a snp
        check_one_exists, check_zero_exists = False, False
        for happ in self.haplotypes_list:
            if self.stype == 'hardsweep' or self.stype == 'softsweep':
                # change non-unity characters
                for i, element in enumerate(happ[:, self.spos_i]):
                    if element[0][0] != '0':
                        happ[i, self.spos_i] = '1'
                # check sweep site is a snp
                if '1' in happ[:, self.spos_i]:
                    check_one_exists = True
                if '0' in happ[:, self.spos_i]:
                    check_zero_exists = True
        # check that sweep site is a snp in all samples
        if check_one_exists == False and check_zero_exists == True:
            self.check = False
            self.problems.add('- snp was lost in all the samples')
        if check_zero_exists == False and check_one_exists == True:
            self.check = False
            self.problems.add('- snp was fixed in all the samples')
        # convert haplotypes from chars to ints
        if self.check:
            self.haplotypes_list = [happ.astype(int) for happ in self.haplotypes_list]
        # # remove spos_i column from haplotypes, posfile, genfile, recmap, varid,
        # if self.remove_spos:
        #     self.varid_list = np.delete(self.varid_list, self.spos_i)
        #     self.physpos_list = np.delete(self.physpos_list, self.spos_i)
        #     self.genpos_list = np.delete(self.genpos_list, self.spos_i)
        #     self.recrates_list = np.delete(self.recrates_list, self.spos_i)
        #     for pop_i, _ in enumerate(self.haplotypes_list):
        #         self.haplotypes_list[pop_i] = np.delete(self.haplotypes_list[pop_i], self.spos_i, axis=1)
        # # # 
        # print problems with simulations
        self.print_problems()
        return self

class HardSweepModel(SweepModel):
    """Class for hard selective sweep models"""
    def __init__(self, params_sims=ParametersSimulations(), key='NA', pop='NA', seed='NA'):
        super(HardSweepModel, self).__init__(params_sims, 'hardsweep', key, pop, seed)

class SoftSweepModel(SweepModel):
    """Class for soft selective sweep models"""
    def __init__(self, params_sims=ParametersSimulations(), key='NA', pop='NA', seed='NA'):
        super(SoftSweepModel, self).__init__(params_sims, 'softsweep', key, pop, seed)
