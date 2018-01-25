#!/usr/bin/env python
# Author: Duncan Tormey
# Email: dut@stowers.org or duncantormey@gmail.com

from __future__ import print_function
from __future__ import division
from itertools import izip
import argparse


def line_to_dict(line, header):
    return dict(zip(header, line.split('\t')))


def key_with_max_val(dictionary):
    values = dictionary.values()
    keys = dictionary.keys()
    return keys[values.index(max(values))]


def ret_allele_1(prof_dict):
    allele_1 = {sorted(prof_dict, key=lambda k: int(prof_dict[k]))[-1]}


    return allele_1

def ret_allele_2(prof_dict, allele_1, fraction):
    nucl = sorted(prof_dict, key=lambda k: prof_dict[k])[-2]
    if prof_dict[nucl] >= fraction *prof_dict[list(allele_1)[0]]:
        allele_2 = {nucl}
    else:
        allele_2 = set()

    return allele_2


def ret_parent_alleles(line, header, fraction):
    prof_dict = line_to_dict(line, header)
    prof_dict = {nuc:int(val) for nuc, val in prof_dict.items() if nuc in ['A','T','C','G'] and int(val) > 0}

    if len(prof_dict.keys()) == 1:
        alleles = {prof_dict.keys()[0]}
    elif len(prof_dict.keys()) >= 2:
        allele_1 = ret_allele_1(prof_dict)
        allele_2 = ret_allele_2(prof_dict, allele_1, fraction)
        alleles = allele_1 | allele_2
    else:
        alleles = set()

    return alleles


def ret_partheno_allele(line, header):
    prof_dict = line_to_dict(line, header)
    prof_dict = {nuc: int(val) for nuc, val in prof_dict.items() if nuc in ['A', 'T', 'C', 'G'] and int(val) > 0}
    if len(prof_dict.keys()) > 0:
        return set(key_with_max_val(prof_dict))
    else:
        return set()


def ret_child_alleles(line, header, fraction, partheno=False):
    if partheno:
        alleles = ret_partheno_allele(line, header)
    else:
         alleles = ret_parent_alleles(line, header, fraction)

    return alleles


def parse_profiles(parent_path, child_path, fraction, partheno=True):
    header = ['chrom', 'pos', 'A', 'C', 'G', 'T', 'cov', 'status']
    num_profs = 0
    num_cont = 0

    child_name = parent_path.split('/')[-1].split('.')[0]
    parent_name = child_path.split('/')[-1].split('.')[0]

    if partheno:
        out_name = './%s_%s_%serr_partheno_prof_diff.tsv' % (child_name, parent_name, str(fraction))
        out_stats_name = './%s_%s_%serr_partheno_prof_diff.stats' % (child_name, parent_name, str(fraction))
    else:
        out_name = './%s_%s_%serr_prof_diff.tsv' % (child_name, parent_name, str(fraction))
        out_stats_name = './%s_%s_%serr_prof_diff.stats' % (child_name, parent_name, str(fraction))

    with open(parent_path) as parent_profs, \
            open(child_path) as child_profs, \
            open(out_name, 'w') as out:

        parent_profs.next()  # skip header
        child_profs.next()  # skip header
        for parent_line, child_line in izip(parent_profs, child_profs):
            num_profs += 1
            parent_alleles = ret_parent_alleles(parent_line, header,fraction)
            child_alleles = ret_child_alleles(child_line, header, fraction, partheno=partheno)

            if child_alleles.issubset(parent_alleles):
                num_cont += 1

            elif len(parent_alleles) == 0:
                num_cont += 1

            elif len(child_alleles) == 0:
                num_cont += 1

            else:
                out.write(parent_line.strip('\n') + '\t' + child_line)

        with open(out_stats_name, 'w') as outstats:
            outstats.write('number_of_explained\t%d\n' % num_cont)
            outstats.write('number_of_sites\t%d\n' % num_profs)
            outstats.write('perc_allele_contribution\t%s\n' % str(float(num_cont) / float(num_profs) * 100))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-parent_path', action='store', dest='parent_path',
                        help='path to file containing parent profiles')

    parser.add_argument('-child_path', action='store', dest='child_path',
                        help='path to file containing child profiles')

    parser.add_argument('-min_perc_reads', action='store', dest='fraction', default=0.6,
                        help='minimum percentage of reads that support an allele for child')

    parser.add_argument('-partheno', action='store_true', dest='partheno', default=False,
                        help='Take only the most abundant allele for the child')

    args = parser.parse_args()

    parse_profiles(args.parent_path, args.child_path, args.fraction, args.partheno)


