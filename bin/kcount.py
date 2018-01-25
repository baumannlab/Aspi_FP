#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com

from __future__ import print_function
import gzip
from collections import Counter, defaultdict
import argparse
import glob
import sys
from functools import wraps
from time import time


def timed(f):
    '''
    This decorator was addapted from: 
    http://www.artima.com/weblogs/viewpost.jsp?thread=240808
   
    It times the execution of whatever function it is applied too.    
    '''
    @wraps(f)
    def wrapper(*args, **kwds):
        start = time()
        result = f(*args, **kwds)
        elapsed = time() - start
        m, s = divmod(elapsed, 60)
        h, m = divmod(m, 60)
        print("%s finished in %d hours, %02d minutes and %02d seconds" % (f.__name__, h, m, s))
        return result

    return wrapper


def ret_reverse_compliment(seq):
    '''
    This function takes a string of DNA sequence
    and retuns the reverse compliment of that sequence

    >>>ret_reverse_compliment('AAAAGGGG')
    'CCCCTTTT'
    '''
    seq = seq[::-1]
    rev_comp = []
    for char in seq:
        if char == 'T':
            rev_comp.append('A')
        elif char == 'A':
            rev_comp.append('T')
        elif char == 'G':
            rev_comp.append('C')
        elif char == 'C':
            rev_comp.append('G')
        else:
            rev_comp.append(char)

    rev_comp = ''.join(rev_comp)
    return rev_comp


def ret_seq_gen(paths):
    for path in paths:
        if '.gz' not in path:
            fq_fh = open(path, 'r')
        else:
            fq_fh = gzip.open(path, 'r')

        for i, line in enumerate(fq_fh):
            if i % 4 == 1:
                yield line.strip()


def ret_kmer_gen(read, ksize):
    for i in xrange(len(read) - ksize + 1):
        yield read[i:i+ksize]


@timed
def count_kmers(paths, ksize):
    kcounter = defaultdict(int)
    for read in ret_seq_gen(paths):
        for kmer in ret_kmer_gen(read, ksize):
            kcounter[kmer] += 1

    return kcounter


@timed
def collapse_strands(kcounter):
    '''
    This function takes a dictionary as input, where the keys are kmers
    and the values are occurances. It iterates through the keys of the
    dictionary and checks if the reverse compliment of the keys are in the
    dictionary, if it is the counts are combined and the reverse compliment
    key is deleted.

    >>> collapse_strands({'ATG':5, 'GGT':5, 'TTA':3, 'TAA':2})
    collapse_strands took 5.48362731934e-05 seconds to finish
    {'ATG': 5, 'TTA': 5, 'GGT': 5}
    '''
    for kmer in kcounter.keys():
        if kmer in kcounter:
            rev_kmer = ret_reverse_compliment(kmer)
            if rev_kmer in kcounter:
                kcounter[kmer] += kcounter[rev_kmer]
                del kcounter[rev_kmer]

    return kcounter


@timed
def ret_multiplicity_data(counter):
    '''
    Creates a counter dictionary from the values of a counter dictionary
    
    >>> ret_multiplicity_data({'ATG':5, 'GGT':5, 'TTA':3, 'TAA':2})
    Counter({5: 2, 2: 1, 3: 1})
    '''
    multi_data = Counter(counter.values())

    return multi_data


@timed
def write_multi_hist(multi_data, outfile):
    '''
    sorts the output of ret_multiplicity_data and writes it to a tsv file
    '''
    fo = open(outfile, 'w')
    multi_data = sorted(multi_data.items(), key=lambda x: x[0])
    fo.write('multiplicity\tnumber_of_kmers\n')
    for m, o in multi_data:
        fo.write('%s\t%s\n' % (str(m), str(o)))
    fo.close()


@timed
def run(args):
    print(args)
    fq_paths = []
    if args.fastq_paths:
        fq_paths.extend(args.fastq_paths)
    if args.wild_card:
        wild_paths = glob.glob(args.wild_card)
        if len(wild_paths) > 0:
            fq_paths.extend(wild_paths)
        else:
            print('wild_card returned no fastq paths')

    if len(fq_paths) == 0:
        print('User must supply at least one fastq file path')
        sys.exit(1)

    kcounter = count_kmers(fq_paths, args.ksize)

    if args.strand:
        print('Strand normalizing...')
        kcounter = collapse_strands(kcounter)

    multi_data = ret_multiplicity_data(kcounter)
    write_multi_hist(multi_data, args.histo)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-ksize', action='store', dest='ksize', type=int,
                        help='size of kmer to use', default=21)

    parser.add_argument('-fastq_paths', action='append', dest='fastq_paths',
                        help='path to fastq file, or paths to fastq files seperated by spaces',
                        default=None)

    parser.add_argument('-wild_card', action='store', dest='wild_card',
                        help='path with wild card used to get list of fastqs',
                        type=str, default=None)
    
    parser.add_argument('-histo', action='store', dest='histo',
                        default='./kmer_histogram.txt',
                        help='name of output histogram tsv file')

    strand_parser = parser.add_mutually_exclusive_group(required=False)

    strand_parser.add_argument('--strand_normalize', dest='strand', action='store_true')

    strand_parser.add_argument('--no_strand_normalize', dest='strand', action='store_false')

    parser.set_defaults(strand=True)
    
    args = parser.parse_args()

    run(args)
