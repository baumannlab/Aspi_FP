#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com
##################################################
# This script takes the tsv file output by      #
# pysam_profiler.py and counts the occurence of #
# each unique profile, ignoring scaffold and    #
# position.                                     #
#################################################


from __future__ import print_function
from __future__ import division
import sys
import os
from collections import Counter


def count_profiles(prof_path):
    outfile = './' + os.path.basename(prof_path).replace('.prof', '.prof_counts')
    counts = Counter()
    with open(prof_path, 'r') as fh:
        fh.next() # skip header
        for line in fh:
            line = line.strip().split()
            prof = tuple(line[2:])
            counts[prof] += 1

    return counts, outfile


def write_profiles(counts, outfile):
    with open(outfile, 'w') as fho:
        fho.write('occurence\tA\tC\tG\tT\tcov\tstatus\n')
        for prof, count in counts.items():
            fho.write('%s\t%s\n' % (str(count), '\t'.join(map(str, prof))))
    

if __name__ == "__main__":
    if len(sys.argv) == 2:
        print(sys.argv)
        counts, outfile = count_profiles(sys.argv[1])
        write_profiles(counts, outfile)
    else:
        print('usage: count_profiles /path/to/profile ')
