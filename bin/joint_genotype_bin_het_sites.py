#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com

from __future__ import print_function
from collections import defaultdict
import numpy as np
import multiprocessing as mp

def get_genotype(splitline, position):
    genotype = splitline[position].split(':')[0].split('/')
    return genotype


def get_chrom(splitline):
    chrom = splitline[0]
    return chrom


def ret_DP_index(splitline):
    index = splitline[8].split(':').index('DP')
    return index


def ret_AD_index(splitline):
    index = splitline[8].split(':').index('AD')
    return index


def ret_coverage(splitline, position):
    try:
        ad_index = ret_AD_index(splitline)
        coverage = int(
            sum([int(x) for x in splitline[position].split(':')[ad_index].split(',')]))
    except ValueError:
        try:
            dp_index = ret_DP_index(splitline)
            coverage = int(splitline[position].split(':')[dp_index])
        except ValueError:
            coverage = 'NoField'
    return coverage


def count_coverage(coverage_count, splitline, position):
    coverage = ret_coverage(splitline, position)

    coverage_count.append(coverage)

    return coverage_count


def is_het(genotype):
    return genotype[0] != genotype[1]


def count_het(animal_count, splitline, position):
    if is_het(get_genotype(splitline, position)):
        animal_count += 1
    else:
        animal_count += 0

    return animal_count


def is_feature(splitline, position):
    feature = eval(splitline[position])
    return feature


def count_feature(feature_count, splitline, position):
    if is_feature(splitline, position):
        feature_count += 1
    else:
        feature_count += 0
    return feature_count


def window_parse(path, animal_names, feature_names, window_size=150):
    animal_counts = defaultdict(int)
    feature_counts = defaultdict(int)
    coverage_dict = defaultdict(list)
    fh = open(path, "r")
    out = path.replace('.vcf', '.%s_feature_window.tsv' % str(window_size))
    fo = open(out, 'w')
    i = 0  # window counter
    k = 0  # chrom checker
    j = 0  # position counter
    temp_chrom = None

    for line in fh:
        if '#CHROM' in line:

            header = line.strip('\n').split()

            animal_positions = [header.index(animal)
                                for animal in animal_names]
            animal_headers = []

            for animal in animal_names:
                animal_headers.extend(
                    [animal + '_het_count', animal + '_het_freq', animal + '_avg_cov', animal + '_std_cov'])

            feature_positions = [header.index(feature)
                                 for feature in feature_names]
            feature_headers = []

            for feature in feature_names:
                feature_headers.extend([feature + '_count', feature + '_perc'])

            header = ['chrom', 'window_start', 'window_stop'] + \
                animal_headers + feature_headers

            print(header)

            fo.write('\t'.join(header) + '\n')

        if line[0] != '#':
            splitline = line.strip('\n').split()
            j += 1  # at a new position

            if k == 0:  # do we know what chromosome we are on?
                chrom = get_chrom(splitline)  # get the chromosome
                k = 1  # we now know what chromosome we are on
                if temp_chrom != chrom:  # is it the chromosome we were on before?
                    j = 0  # no so start the position counter over

            # have we filled the window and is this line part of the current
            # windows chromosome?
            if i < window_size and chrom == get_chrom(splitline):
                i += 1  # window counter increases
                for position in animal_positions:  # look at genotypes for each sample
                    # count if the genotype is heterozygous
                    animal_counts[position] = count_het(
                        animal_counts[position], splitline, position)

                    coverage_dict[position] = count_coverage(
                        coverage_dict[position], splitline, position)

                for position in feature_positions:
                    feature_counts[position] = count_feature(
                        feature_counts[position], splitline, position)

            else:  # Oh snap the window is full, or we switched chromosomes
                # store the begining of our data line
                outline = [chrom, j - i, j]
                temp_chrom = chrom  # keep what chromosome we were on

                for position in animal_positions:  # look at het counts for each sample
                    # make some frequencies and add them to our out line
                    outline.extend([animal_counts[position],
                                    float(animal_counts[position]) / float(i),
                                    np.mean(coverage_dict[position]),
                                    np.std(coverage_dict[position])])
                    animal_counts[position] = 0  # reset het counts
                    coverage_dict[position] = []  # rest coverage_dict

                for position in feature_positions:
                    outline.extend([feature_counts[position],
                                    float(feature_counts[position]) / float(i)])

                    feature_counts[position] = 0

                outstrings = [str(x)
                              for x in outline]  # string format out line
                fo.write("\t".join(outstrings) + "\n")  # write the out line
                # we no longer no what chromosome is next(could be same chrom
                # next window, or new chrom next window)
                k = 0
                i = 0  # window counter gets reset.
    fh.close()
    fo.close()
    return out

if __name__ == "__main__":
    paths = ['/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/final_variant_calling/joint_genotypes/jg_A_tigris8450_Atig_122_Atig003_Atig001.no_call_removed.filt.recode.annotated.norep.vcf',
             '/home/dut/projects/tigris/heterozygosity/dwn_sample_atig_122/final_variant_calling/joint_genotypes/jg_A_tigris8450_Atig_122_Atig003_Atig001.no_call_removed.filt.recode.annotated.vcf']

    animals = ["A_tigris8450", "Atig001", "Atig003", "Atig_122"]
    features = ["intron", "five_prime_UTR", "exon", "mRNA", "CDS", "three_prime_UTR", "gene", "dispersed_repeat"]
    sizes = [10000, 1000, 150]
    cores = 3
    pool = mp.Pool(processes=cores)
    results = [pool.apply_async(window_parse, args=(paths[1], animals, features, size)) for size in sizes]
    output = [p.get() for p in results]
