#--------------------------------------------------------#
# Stephen Rong, Ramachandran Lab, Brown University       #
# Created: 2016-12-06, Updated: 2017-02-12               # 
# Descrption: Preprocess 1000 Genomes Phase 1 haplotypes # 
#--------------------------------------------------------#
#!/usr/bin/env python

from __future__ import division
import numpy as np
import pandas as pd
from pygg import *
from subprocess import Popen, PIPE
from remove_relations import remove_cryptic_relations

def get_unrelateds_temp(individuals_in, relations_in, relations_retained, relations_removed, relations_dot, relations_fig):
    """Returns dict of populations and their individuals with cryptic relations removed."""
    # get sampled individuals
    individuals_list = set()
    with open(individuals_in, 'rb') as f:
        next(f)
        for i, line in enumerate(f):
            ind_temp = line.split()
            individuals_list.add(ind_temp[0])
    # get list of cryptic relations
    relations_list = set()
    with open(relations_in, 'rb') as f:
        next(f)
        for i, line in enumerate(f):
            relations_temp = line.split()
            relid, relpop = relations_temp[1], relations_temp[4]
            relrels = [i.split(':')[0] for i in relations_temp[7].split(';')]
            #    only save cryptic relations_in if both individuals_in present
            for j in relrels:
                if relid in individuals_list and j in individuals_list:
                    relations_list.add((relid, j))
    # get set of removed individauls
    removed_list = remove_cryptic_relations(relations_list, relations_dot, relations_fig)
    # get unrelated individuals by population
    individuals_temp = dict()
    with open(individuals_in, 'rb') as f:
        with open(relations_retained, 'w') as g:
            with open(relations_removed, 'w') as h:
                next(f)
                # iterate individuals
                for i, line in enumerate(f):
                    ind_temp = line.split()
                    indid, indpop, indreg, indsex = ind_temp[0], ind_temp[1], ind_temp[2], ind_temp[3]
                    # retained individuals
                    if indid not in removed_list:
                        if indpop not in individuals_temp:
                            individuals_temp[indpop] = set()
                        individuals_temp[indpop].add(i)
                        g.write('{0} {1} {2} {3} {4}\n'.format(i, indid, indpop, indreg, indsex))
                    # removed individuals
                    else:
                        h.write('{0} {1} {2} {3} {4}\n'.format(i, indid, indpop, indreg, indsex))
    # convert set into sorted list of individuals
    for pop in individuals_temp:
        individuals_temp[pop] = sorted(list(individuals_temp[pop]))
    return individuals_temp

def get_genetic_positions(chrm, legend_temp, genetic_map_in, genpos_fig, recmap_fig, cytoband_in, sgp, srm):
    """Return dict of variant indices with interpolated genetic positions, 
        and plot genetic position and recombination rate vs physical position."""
    # get input files
    genetic_temp = pd.read_csv(genetic_map_in, sep = ' ', header = 0, names = ['physpos', 'recrate', 'genpos'])
    cytoband_temp = pd.read_csv(cytoband_in, sep = '\t', header = 0, names = ['chrm', 'lower', 'upper', 'locus', 'band'])
    centromere_temp = cytoband_temp[(cytoband_temp.band == 'acen') & (cytoband_temp.chrm == ('chr'+chrm))]
    # linear interpolation
    x_physpos = np.array(genetic_temp.physpos)
    y_genpos = np.array(genetic_temp.genpos)
    z_recrate = np.array(genetic_temp.recrate)
    # get physical position
    x = np.zeros(len(legend_temp))
    legend_sort = sorted(legend_temp)
    for i, var in enumerate(legend_sort):
        x[i] = legend_temp[var][1]
    # get genetic position
    y = np.interp(x, x_physpos, y_genpos)
    # plot genetic position
    p = ggplot(dict(x=x[::sgp], y=y[::sgp]), aes(x='x', y='y')) + geom_line()
    for i, row in centromere_temp.iterrows():
        p = p + geom_vline(xintercept = row[1], colour = 2) + geom_vline(xintercept = row[2], colour = 2)
    ggsave(genpos_fig, p)
    # get recombination rate
    z = np.interp(x[::srm], x_physpos[::srm], z_recrate[::srm])
    # plot recombination map
    q = ggplot(dict(x=x[::srm], z=z), aes(x='x', y='z')) + geom_line()
    for i, row in centromere_temp.iterrows():
        q = q + geom_vline(xintercept = row[1], colour = 2) + geom_vline(xintercept = row[2], colour = 2)
    ggsave(recmap_fig, q)
    # store genetic position
    for i, var in enumerate(legend_sort):
        genpos = y[i]
        recrate = z[i]
        legend_temp[var].insert(2, genpos)
        legend_temp[var].insert(3, recrate)
    return legend_temp

def get_legend_temp(chrm, legend_gz, genetic_map_in, genpos_fig, recmap_fig, cytoband_in, sgp, srm):
    """Return dict of variant indices of SNPs only, checking for concurrent SNPs, 
        and filling in missing varids based on plink2 --set-mising-var-ids."""
    legend_temp = dict()
    content = Popen(args = 'zcat < {0}'.format(legend_gz), stdout=PIPE, shell=True).communicate()[0].strip().split('\n')
    physpos_last = False
    # shift by one for header
    for i, line in enumerate(content[1:]):
        var_temp = line.split()
        varid, physpos, vartype = var_temp[0], var_temp[1], var_temp[4]
        # filter out non-SNPs
        if vartype == 'SNP':
            # check for concurrent SNPs
            if int(physpos_last) == int(physpos):
                assert False, 'Concurrent SNPs found at position {0}'.format(physpos_last)
            physpos_last = physpos
            # fill in missing variant IDs
            if varid == '.':
                varid = 'chr{0}:{1}[b37]{2},{3}'.format(chrm, physpos, var_temp[2], var_temp[3])
            # set placeholder without genpos
            legend_temp[i] = [varid, physpos]  # , vartype]
    legend_temp = get_genetic_positions(chrm, legend_temp, genetic_map_in, genpos_fig, recmap_fig, cytoband_in, sgp, srm)
    return legend_temp

def out_hapfiles(chrm, individuals_temp, legend_temp, haplotypes_gz, hap_files):
    """Generates hap_files for each chromosome in each population being considered, 
        while also filtering out monomorphic variants w.r.t. populations being considered."""
    content = Popen(args = 'zcat < {0}'.format(haplotypes_gz), stdout=PIPE, shell=True).communicate()[0].strip().split('\n')
    # open hapfile for population
    hap_f = dict()
    # get haplotype indices for population
    hap_temp = dict()
    for pop in hap_files:
        hap_f[pop] = open(hap_files[pop], 'w')
        hap_temp[pop] = np.asarray(individuals_temp[pop])
        hap_temp[pop] = np.concatenate([2*hap_temp[pop], 2*hap_temp[pop]+1])
    hap_temp['NA'] = np.concatenate([hap_temp[pop] for pop in hap_files])
    # sort indices so they appear in order
    for pop in hap_temp:
        hap_temp[pop] = np.sort(hap_temp[pop])
    # go through each variant/line in haplotypes_gz
    hap_total = np.size(hap_temp['NA'])
    for i, line in enumerate(content):
        # check if variant in legend
        if i in legend_temp:
            np_line = np.fromstring(line, dtype = np.int8, sep = ' ')
            np_line_sum = np.sum(np_line[hap_temp['NA']])
            # check if monomorphic with respect to included populations            
            if np_line_sum > 0 and np_line_sum < hap_total:
                for pop in hap_files:
                    np_line_pop = np_line[hap_temp[pop]]
                    hap_f[pop].write(' '.join(np_line_pop.astype(np.str)) + '\n')
            # remove variant from analysis if monomorphic in populations
            else:
                del legend_temp[i]
    # close hapfile for population
    for pop in hap_files:
        hap_f[pop].close()
    return legend_temp

def out_posfiles(chrm, legend_temp, pos_file, pos_recfile):
    """Output .pos file for selscan or hapbin."""
    with open(pos_file, 'w') as f:
        with open(pos_recfile, 'w') as g:
            for i in sorted(legend_temp.keys()):
                varid = legend_temp[i][0]
                physpos = legend_temp[i][1]
                genpos = legend_temp[i][2]
                recrate = legend_temp[i][3]
                g.write('{0}\n'.format(recrate))
                f.write('{0} {1} {2} {3}\n'.format(chrm, varid, genpos, physpos))

if __name__ == '__main__':
    # populations
    pop_list = ['YRI', 'CEU', 'CHB']  # # # 
    # directories
    prefix_raw = '../../data/applications/1000GP_Phase1_raw/'
    prefix_post = '../../data/applications/1000GP_Phase1_post/'
    prefix_misc = '../../data/applications/1000GP_Phase1_misc/'
    # input files
    individuals_in = prefix_raw + 'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample'
    relations_in = prefix_misc + '20120220_cryptic_relations_clusters.txt'
    cytoband_in = prefix_misc + 'cytoband.txt'
    # output files
    relations_retained = prefix_post + 'retained_individuals.txt'
    relations_removed = prefix_post + 'removed_individuals.txt'
    relations_dot = prefix_post + 'cryptic_relations_graph.dot'
    relations_fig = prefix_post + 'cryptic_relations_graph.pdf'
    # iterate  chromosomes
    for i in range(1, 23)[::-1]:
        chrm = str(i)
        print 'Chromosome {0}...'.format(chrm)
        # input files per chromosome
        legend_gz = prefix_raw + 'ALL.chr{0}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz'.format(chrm)
        haplotypes_gz = prefix_raw + 'ALL.chr{0}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz'.format(chrm)
        genetic_map_in = prefix_raw + 'genetic_map_chr{0}_combined_b37.txt'.format(chrm)
        # output files per chromsome
        genpos_fig = prefix_post + 'genpos_chr{0}.pdf'.format(chrm)
        recmap_fig = prefix_post + 'recmap_chr{0}.pdf'.format(chrm)
        pos_file = prefix_post + '1000GPP1_chr{0}.pos'.format(chrm)
        pos_recfile = prefix_post + '1000GPP1_chr{0}_rec'.format(chrm)
        # subsample genpos and recmap figures
        sgp, srm = 1, 1  # no subsample: 1, 1
        # store haplotype file locations in dict
        hap_files = dict()
        for pop in pop_list:
            hap_files[pop] = prefix_post + '1000GPP1_chr{0}_pop{1}.hap'.format(chrm, pop)
        # output hapfiles and posfile per chromosome
        #     remove cryptic relations
        individuals_temp = get_unrelateds_temp(individuals_in, relations_in, relations_retained, relations_removed, relations_dot, relations_fig)
        #     remove non-SNP variants, fill in variant ids, interpolate genetic positions
        legend_temp = get_legend_temp(chrm, legend_gz, genetic_map_in, genpos_fig, recmap_fig, cytoband_in, sgp, srm)
        #     remove monomorphic variants among included pops, output hapfiles for each pop
        legend_temp = out_hapfiles(chrm, individuals_temp, legend_temp, haplotypes_gz, hap_files)
        #     output final posfile, used with hapfiles for selscan or hapbin
        out_posfiles(chrm, legend_temp, pos_file, pos_recfile)
