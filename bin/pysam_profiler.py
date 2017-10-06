#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com
#################################################
# This script takes a bam file and a reference  #
# genome file as input. It uses the pysam       #
# python package to pull out the per position   #
# nucleotide coverage, and determes the statust #
# of the position(Transition, Transversion,     #
# Homozygous or Unknown)                        #
#################################################

from __future__ import print_function
from __future__ import division
import pysam
import pysamstats
import sys
import os

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

    prof.extend([sum(prof), status])

    return  prof

def write_profiles(bam_path, ref_path):
    outfile = './' + os.path.basename(bam_path).replace('.bam', '.prof')
    bam = pysam.AlignmentFile(bam_path)
    with open(outfile, 'w') as fh:
        fh.write('chrom\tpos\tA\tC\tG\tT\tcov\tstatus\n')
        for rec in pysamstats.stat_variation(bam, ref_path, pad=True, max_depth=1000000):
            prof = [rec['A'], rec['C'], rec['G'], rec['T']]
            prof = get_prof_status(prof)
            prof[:0] = [rec['chrom'], rec['pos']]
            fh.write('%s\n' % '\t'.join(map(str,prof)))


if __name__ == "__main__":
    if len(sys.argv) == 3:
        print(sys.argv)
        write_profiles(sys.argv[1], sys.argv[2])
    else:
        print('usage: pysam_profiler.py /path/to/alignment.bam /path/to/reference.fa')
