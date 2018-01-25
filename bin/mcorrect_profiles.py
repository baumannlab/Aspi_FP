#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com
##################################################
# This function sets any allele at a heterozygous#
# or unknown site in the parthenogenetic animal  #
# to zero if it is not present in the mother     #
# genotype at that position.                     #
from __future__ import print_function
from __future__ import division
import pandas as pd
import multiprocessing as mp
import numpy as np


def assign_string_type(string):
    if not string:
        value = None
        
    if '.' in string:
        try:
            value = float(string)
        except ValueError:
            value = string
    else:
        try:
            value = int(string)
        except ValueError:
            value = string

    return value


def get_prof_status(prof, bases=['A', 'C', 'G', 'T']):
    bases = [bases[i] for i, n in enumerate(prof) if n > 0]
    if len(bases) == 2:
        if bases == ['A', 'G'] or bases == ['C', 'T']:
            status = 'Transition'
        else:
            status = 'Transversion'
    elif len(bases) == 1:
            status = 'Homozygous'
    else:
        status = 'Unknown'


    return  status

def ret_pos_freq_dict(path):
    pos_freqs = {}
    with open(path,'r') as fh:
        fh.next()
        for line in fh:
            line=line.strip('\n').split()
            pos_freqs[(line[0], int(line[1]), 'A')] = int(line[2])
            pos_freqs[(line[0], int(line[1]), 'C')] = int(line[3])
            pos_freqs[(line[0], int(line[1]), 'G')] = int(line[4])
            pos_freqs[(line[0], int(line[1]), 'T')] = int(line[5])

    return pos_freqs


def ret_dict_line(header,line):
    line = [assign_string_type(field) for field in line.strip('\n').split()]
    return dict(zip(header,line))


def correct_prof_dict(prof_dict, pos_freqs):
    nucs = ['A','C','G','T']
    mnucs = ['mA','mC','mG','mT']
    for n,m in zip(nucs, mnucs):
        if prof_dict[n] != 0 and pos_freqs[(prof_dict['chrom'], prof_dict['pos'], n)] == 0:
            prof_dict[m] = 0
        else:
            prof_dict[m] = prof_dict[n]

            
    prof_dict['mcov'] = sum((prof_dict[m] for m in mnucs))
    prof_dict['mstatus'] = get_prof_status([prof_dict[m] for m in mnucs])

    return prof_dict


def write_corrected_daughter(daughter_path, mother_path):
    pos_freqs = ret_pos_freq_dict(mother_path)
    corrected_path = daughter_path.replace('.prof','.mcor.prof')
    with open(daughter_path, 'r') as fh, open(corrected_path, 'w') as fho:
        header = fh.readline().strip('\n').split()
        new_header = sorted(header[0:2] + ['m%s' % h for h in header[2:]] + header[2:])

        fho.write('%s\n' % '\t'.join(map(str, new_header)))
        
        for line in fh:
            prof_dict = ret_dict_line(header, line)
            prof_dict = correct_prof_dict(prof_dict, pos_freqs)
            outline = '%s\n' % '\t'.join([str(prof_dict[h]) for h in new_header])
            fho.write(outline)

if __name__ == '__main__':
    mother = '/home/dut/projects/tigris/heterozygosity/pysam/Atig_122.merged.dedup.realigned.prof'
    daughter = '/home/dut/projects/tigris/heterozygosity/pysam/A_tigris8450.merged.dedup.realigned.het.prof'
    write_corrected_daughter(daughter,mother)
