#!/usr/bin/env python
# Author: Duncan Tormey
# Email: dut@stowers.org or duncantormey@gmail.com
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
import argparse


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


def prof_generator(prof_path):
    with open(prof_path, 'r') as fh:
        header = fh.next().strip('\n').split('\t')
        for line in fh:
            prof = [assign_string_type(l) for l in line.strip('\n').split('\t')]
            prof_dict = dict(zip(header, prof))
            yield prof_dict


def write_out_line(prof_dict, sum_coverage, window_size, fh):
    start = prof_dict['pos'] - window_size
    if start < 0:
        start += 1
    stop = prof_dict['pos']
    outline = '%s\t%s\t%s\t%s\t%s\n' % (prof_dict['chrom'],
                                        start, stop, sum_coverage, window_size)
    fh.write(outline)


def write_coverage_windows(prof_path, window=10000):
    out_file_name = prof_path.replace('.prof', '.%skb.cov.windows.tsv' % str(window / 1000))
    file_out = open(out_file_name, 'w')
    print(out_file_name)
    prof_gen = prof_generator(prof_path)
    first_prof = next(prof_gen)
    window_size = 0
    sum_coverage = first_prof['cov']
    current_scaffold = first_prof['chrom']
    for prof_dict in prof_generator(prof_path):
        window_size += 1
        if window_size <= window and current_scaffold == prof_dict['chrom']:
            sum_coverage += prof_dict['cov']

        elif current_scaffold == prof_dict['chrom']:
            write_out_line(prev_prof_dict, prev_sum_coverage, prev_window_size, file_out)
            sum_coverage = prof_dict['cov']
            window_size = 0
            current_scaffold = prof_dict['chrom']
        else:
            write_out_line(prev_prof_dict, prev_sum_coverage, prev_window_size, file_out)
            sum_coverage = prof_dict['cov']
            window_size = 0
            current_scaffold = prof_dict['chrom']

        prev_sum_coverage = sum_coverage
        prev_prof_dict = prof_dict
        prev_window_size = window_size
    file_out.close()

    return out_file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-prof', action='store', dest='prof_path',
                        help='input profile')

    parser.add_argument('-window_size', action='store', dest='window_size', default=10000, type=int,
                        help='designate the window size (default 10000)')
    args = parser.parse_args()

    write_coverage_windows(args.prof_path, window=args.window_size)
