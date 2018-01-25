#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com

from __future__ import print_function
from collections import defaultdict
import multiprocessing as mp
import numpy as np

def get_genotype(splitline, position):
    genotype = splitline[position].split(':')[0].split('/')
    return genotype


def get_DP_index(splitline):
    index = splitline[8].split(':').index('DP')
    return index


def get_AD_index(splitline):
    index = splitline[8].split(':').index('AD')
    return index


def ret_coverage(splitline, position):
    try:
        ad_index = get_AD_index(splitline)
        coverage = int(
            sum([int(x) for x in splitline[position].split(':')[ad_index].split(',')]))
    except ValueError:
        try:
            dp_index = get_DP_index(splitline)
            coverage = int(splitline[position].split(':')[dp_index])
        except ValueError:
            coverage = 'NoField'
    return coverage


def is_het(genotype):
    return genotype[0] != genotype[1]


def window_parse(path, animal_names):
    fh = open(path, "r")
    out = path.replace('.vcf', '.cov_ratio.tsv')
    fo = open(out, 'w')
    list_structure = np.array([[0, 0]]*len(animal_names))
    cov_dict =  {}
    l = 0
    for line in fh:
        l += 1
        if '#CHROM' in line:
            line = line.strip('\n').split()
            animal_positions = [line.index(a) for a in animal_names]
            print(zip(animal_names, animal_positions))
        elif line[0] != '#':
            splitline = line.strip('\n').split()
            for i in range(len(animal_positions)):
                coverage = ret_coverage(splitline, animal_positions[i])
                genotype = get_genotype(splitline, animal_positions[i])
                if coverage not in cov_dict:
                    structure_copy = list_structure.copy()
                    cov_dict[coverage] = structure_copy
                    
                if is_het(genotype):
                    
                    cov_dict[coverage][i][0] += 1
                else:
                    cov_dict[coverage][i][1] += 1
    header = ['coverage']
    print(l)
    for animal in animal_names:
        header.extend(['%s_het' % animal, '%s_hom' % animal])
    print(header)
    header = '\t'.join(header) + '\n'
    fo.write(header)
    for key in cov_dict:
        flat = [str(item) for sublist in cov_dict[key] for item in sublist]
        line = str(key) + '\t' + '\t'.join(flat) + '\n'
        fo.write(line)
    fh.close()
    fo.close()

    return out

if __name__ == "__main__":
    paths = ['/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/final_variant_calling/joint_genotypes/jg_A_tigris8450_Atig_122_Atig003_Atig001.no_call_removed.filt.recode.annotated.norep.vcf',
             '/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/final_variant_calling/joint_genotypes/jg_A_tigris8450_Atig_122_Atig003_Atig001.no_call_removed.filt.recode.annotated.vcf']
    animals = ["A_tigris8450", "Atig001", "Atig003", "Atig_122"]
    cores = 2
    pool = mp.Pool(processes=cores)
    results = [pool.apply_async(window_parse, args=(path, animals)) for path in paths]
    output = [p.get() for p in results]
