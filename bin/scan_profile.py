#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com
##################################################
# This script takes the  tsv output from         #
# pysam_profiler.py and applies a 10Kb sliding   #
# window accross each scaffod, counting sites as #
# heterozygous if the pass the conditions of the #
# function is_het. This script is currently hard #
# coded for my data.                             #
##################################################

from __future__ import print_function
from __future__ import division
import pandas as pd
import multiprocessing as mp
import numpy as np

def apply_df(df, func, *args):
    return df.apply(lambda x: func(x, *args), axis=1)


def apply_by_multiprocessing(df, func, workers, *args):
    pool = mp.Pool(processes=workers)
    result = [pool.apply_async(apply_df, args = (d, func) + args) for d in np.array_split(df, workers)]
    output = [p.get() for p in result]
    pool.close()
    return pd.concat(output)



def check_equal_het(row):
    counts = [row['A'],row['T'],row['C'],row['G']]
    non_zero_counts = [n for n in counts if n > 0]
    if len(non_zero_counts)==2 and len(set(non_zero_counts))==1:
        return True
    else:
        return False


def is_het(row, avg_cov):
    if row['cov'] > avg_cov*2 and row['cov']<=8:
        return False
    elif check_equal_het(row):
        return True
    else:
        return False


def load_het_profiles(path, cpus, avg_cov):
    het_profs = pd.read_csv(path,sep='\t',header=0)
    het_profs['is_het'] = apply_by_multiprocessing(het_profs, is_het, cpus, avg_cov) #this broke with all sites
    #het_profs['is_het'] = het_profs.apply(lambda row: is_het(row, avg_cov))
    het_profs[het_profs['is_het']==True].to_csv(path + '.is_het8.2avg',sep = '\t', index=False)
    
    return het_profs


def return_window_df(het_profs, scaffold_sizes, window=10000):
    windows = []
    append = windows.append
    for scaffold in het_profs.chrom.unique():
        scaffold_data = het_profs[het_profs.chrom == scaffold]
        scaffold_size = scaffold_sizes[scaffold_sizes.scaffold == scaffold]['scaffold_size'].values[0]
        for i in xrange(0,scaffold_size, window):
            if i+window <= scaffold_size:
                num_het = scaffold_data[(scaffold_data['pos']>=i)&(scaffold_data['pos']<i+window)]['is_het'].sum()
                size = window
            else:
                num_het = scaffold_data[scaffold_data['pos']>i]['is_het'].sum()
                size = len(scaffold_data[scaffold_data['pos']>i])

            append((scaffold, i, num_het, size))

    window_df = pd.DataFrame(windows, columns = ['chrom','window_start', 'het_sites','window_size'])
    
    return window_df
    
def write_window_df(het_prof_path, avg_cov, scaffold_sizes_path, cpus, window=10000):
    scaffold_sizes = pd.read_csv(scaffold_sizes_path,sep='\t', names = ['scaffold','scaffold_size'])
    het_profs = load_het_profiles(het_prof_path, cpus, avg_cov)
    window_df = return_window_df(het_profs, scaffold_sizes, window)
    window_df.to_csv(het_prof_path+'.%skbwindows8.tsv' % str(window/1000), sep='\t',index=False)
    return None

if __name__ == '__main__':
    scaffold_sizes_path = '/home/dut/projects/tigris/genome_annotation/fasta/scaffold_sizes.clean.tsv'

    paths =['/n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/pysam/Atig001.merged.dedup.realigned.prof.not_hom',
            '/n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/pysam/Atig003.merged.dedup.realigned.prof.not_hom',
            '/n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/pysam/Atig_122.merged.dedup.realigned.prof.not_hom',
            '/n/projects/dut/a_marmorata/parthenogen_heterozygosity/data/pysam/A_tigris8450.merged.dedup.realigned.prof.not_hom']

    avg_covs = [16, 18,19,18]
    for path, avg_cov in zip(paths, avg_covs):
        write_window_df(path, avg_cov, scaffold_sizes_path, 8, 10000)

    
