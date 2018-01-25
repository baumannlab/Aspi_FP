#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com

from __future__ import print_function
import datetime
import os
import gc
import sys
import cPickle as pickle
import argparse

gc.collect()
del gc.garbage[:]


def read_list_file(path):
    l = list(set(open(path.strip('\n'), 'r').read().strip('\n').splitlines()))

    return l


def ret_genotype(splitline, position):
    genotype = splitline[position].split(':')[0]
    return genotype


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


def is_het(genotype):
    return genotype.split('/')[0] != genotype.split('/')[1]


def ret_intervals_dict(intervals_path, scaffold_column, feature_column,
                       start_column, end_column, intervals_dict={}):
    i = 0
    preview = range(10)
    intervals_file = open(intervals_path, 'r').readlines()
    print('Reading intervals file %s' % intervals_path)
    for line in intervals_file:
        if line[0] == '#':
            continue
        line = line.strip('\n').split()
        scaffold_name = line[scaffold_column]
        start_coord = int(line[start_column])
        end_coord = int(line[end_column])
        interval = set(range(start_coord, end_coord + 1))
        feature = line[feature_column]

        if (scaffold_name, feature) not in intervals_dict:
            intervals_dict[(scaffold_name, feature)] = interval

        else:
            intervals_dict[(scaffold_name, feature)].update(interval)

        if i in preview:
            print(line)
            print('scaffold_name = %s, start_coord = %s,'
                  ' end_coord = %s, feature = %s' % (scaffold_name,
                                                     start_coord,
                                                     end_coord, feature))
        if i == 10:
            print('...')
        if i % 100000 == 0:
            print(i)
        i += 1
    return intervals_dict


def annotate_vcf_file(intervals_dict, vcf_path, animal_ids):
    out_name = '.'.join([x for x in vcf_path.split('.')[:-1]]) + \
        '.annotated.' + vcf_path.split('.')[-1]
    annot_list = list(set([x[1] for x in intervals_dict.keys()]))
    fh = open(vcf_path, 'r')
    fho = open(out_name, 'w')
    start = datetime.datetime.now()
    print('Annotating vcf file: %s' % vcf_path)
    i = 0
    for line in fh:
        if '#' not in line:
            line = line.strip('\n').split()
            for position in animal_positions:
                line.append(ret_genotype(line, position))
                line.append(ret_coverage(line, position))
            for annot in annot_list:
                if (line[0], annot) in intervals_dict and int(line[1]) in intervals_dict[(line[0], annot)]:
                    line.append('True')

                else:
                    line.append('False')

            fho.write('\t'.join([str(s) for s in line]) + '\n')
        elif '#CHROM' in line:
            line = line.strip('\n').split()
            # add_fields = ['genotype_%s' % x, 'coverage_%s' % x for x in
            # animal_ids]
            add_fields = [s % x for x in animal_ids for s in (
                'genotype_%s', 'coverage_%s')]
            print(add_fields)
            animal_positions = [line.index(a) for a in animal_ids]
            line.extend(add_fields)
            line.extend(annot_list)
            fho.write('\t'.join(line) + '\n')

        else:
            fho.write(line)

        i += 1
        if i % 1000000 == 0:
            print(i)

    print(datetime.datetime.now() - start)
    fh.close()
    fho.close()


def run(args):
    print(args)
    wrkdir = os.getcwd()
    gff3_list = read_list_file(args.gff3_file_list)
    animal_ids = read_list_file(args.animal_id_file_list)
    pickle_name = '%s/%s.pickle' % (wrkdir,
                                    '_'.join([os.path.basename(x).replace('.gff3', '') for x in gff3_list]))
    if os.path.exists(pickle_name):
        intervals_dict = pickle.load(pickle_name, 'rb')
    else:
        intervals_dict = {}
        for gff in gff3_list:
            intervals_dict = ret_intervals_dict(
                intervals_path=gff,
                scaffold_column=args.scaffold_column,
                feature_column=args.feature_column,
                start_column=args.start_column,
                end_column=args.end_column, intervals_dict=intervals_dict)

        print('Pickling intervals_dict....')
        outpickle = open(pickle_name, 'wb')
        pickle.dump(intervals_dict, outpickle)
        print('Done pickling intervals_dict, pickle saved to %s' % pickle_name)

    annotate_vcf_file(intervals_dict, args.vcf_path, animal_ids=animal_ids)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-gff3_file_list', action='store', dest='gff3_file_list',
                        help='path to text file containing all gff3 file paths'
                        'to be used in execution, seperated by newline breaks')

    parser.add_argument('-animal_id_file_list', action='store', dest='animal_id_file_list',
                        help='path to text file containing all animal_ids  in '
                        'order of appearce as they occur in the vcf file'
                        ' to be used in execution, seperated by newline breaks')

    parser.add_argument('-vcf_path', action='store', dest='vcf_path',
                        help='path to vcf file to be annotated')

    parser.add_argument('-start_coord_column', action='store', dest='start_column', default=3,
                        help='position of column containing start coordinate (python counting)')

    parser.add_argument('-end_coord_column', action='store', dest='end_column', default=4,
                        help='position of column containing end coordinate (python counting)')

    parser.add_argument('-fearure_column', action='store', dest='feature_column', default=2,
                        help='position of column containing feature_name (python counting)')

    parser.add_argument('-scaffold_column', action='store', dest='scaffold_column', default=0,
                        help='position of column containing sscaffold or chromosome name (python counting)')
    args = parser.parse_args()
    run(args)
