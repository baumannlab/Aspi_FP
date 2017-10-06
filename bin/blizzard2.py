#!/usr/bin/env python
# Author: Duncan Tormey
# Email: dut@stowers.org or duncantormey@gmail.com

from __future__ import print_function
from __future__ import division
from collections import Counter
import pandas as pd


class MicrosatellitePopulation(object):
    """
    Class for calculating iternal relatedness and homozygosity by loci for a set of individuals
    micosatellite genotypes. Also provides a a way to filter data prior to analysis.
    """

    def __init__(self, data_file_path, split_id=True):

        self.data_file_path = data_file_path
        self.genotypes_df = pd.read_excel(self.data_file_path)
        if split_id:
            self.genotypes_df.columns = ['sample_id', 'size_1', 'size_2', 'height_1', 'height_2']
            self.genotypes_df['sample_name'] = self.genotypes_df['sample_id'].apply(lambda x: x.split('-')[0])
            self.genotypes_df['micro_sat'] = self.genotypes_df['sample_id'].apply(lambda x: x.split('-')[1])
        else:
            self.genotypes_df.columns = ['sample_name', 'micro_sat', 'size_1', 'size_2', 'height_1', 'height_2']

        self.genotypes_df.size_2.fillna(self.genotypes_df.size_1, inplace=True)
        self.genotypes_df.drop_duplicates(inplace=True)
        self.population_sizes = {}
        self.sat_bins = {}
        self.population_intervals = {}
        self.binned_population_sizes = {}
        self.population_size_frequencies = {}


    def remove_microsattelite(self, micro_sat):
        """Removes a specific microsatellite from all individuals in population data"""
        self.genotypes_df = self.genotypes_df[self.genotypes_df.micro_sat != micro_sat]

    def remove_sample(self, sample_name):
        """Removes a individual based on sample name from the population"""
        self.genotypes_df = self.genotypes_df[self.genotypes_df.sample_name != sample_name]
        self.genotypes_df = self.genotypes_df[self.genotypes_df.sample_name != str(sample_name)]

    def get_population_sizes(self):
        """Populates the population sizes dictionary with all alleles across all animals"""
        for sat in self.genotypes_df.micro_sat.unique():
            size_dist = \
                self.genotypes_df[self.genotypes_df.micro_sat == sat].filter(regex='size').stack().reset_index()[
                    0].tolist()
            size_dist = [round(size) for size in size_dist]
            self.population_sizes[sat] = size_dist

    def get_sat_bins(self):
        """Determines the number of bins for each microsatellite, based on 3 nucleotide windows"""
        for sat in self.population_sizes:
            bins = round(float(max(self.population_sizes[sat]) -
                               min(self.population_sizes[sat])) / 3.0)
            if bins <= 0:
                bins = 1
            self.sat_bins[sat] = bins

    def get_population_intervals(self):
        """Determines the intervals for the number of bins for each microsatellite based on
        self.get_sat_bins"""
        for sat in self.population_sizes:
            out, bins = pd.cut(self.population_sizes[sat],
                               bins=self.sat_bins[sat],
                               retbins=True)

            intervals = [
                tuple(
                    float(i)
                    for i in str(o).translate(None, "()[]").split(', ')
                )
                for o in set(out)
                ]
            self.population_intervals[sat] = intervals

    def in_interval(self, size, intervals):
        """This is method that takes a allele and list of intervals and returns the interval in
        which the allele resides"""
        matched = None
        for interval in intervals:
            if interval[0] < size <= interval[1]:
                matched = interval[1]
                break
        return matched

    def get_binned_population_sizes(self):
        """Assigns each of the alleles in population_sizes to an interval"""
        for sat in self.population_sizes:
            self.binned_population_sizes[sat] = [self.in_interval(size, self.population_intervals[sat])
                                                 for size in self.population_sizes[sat]]

    def get_population_size_frequencies(self):
        """Determines the frequency of each allele in binned_population_sizes."""
        self.population_size_frequencies = {
            key: {
                k: float(v) / float(len(val))
                for k, v in Counter(val).items()
                }
            for key, val in self.binned_population_sizes.items()
            }

    def bin_data(self):
        """Bins thea data for each sample"""
        self.get_population_sizes()
        self.get_sat_bins()
        self.get_population_intervals()
        self.get_binned_population_sizes()
        self.get_population_size_frequencies()
        self.genotypes_df['binned_size_1'] = self.genotypes_df.apply(lambda x:
                                                                     self.in_interval(round(x['size_1']),
                                                                                      self.population_intervals[
                                                                                          x['micro_sat']]),
                                                                     axis=1)
        self.genotypes_df['binned_size_2'] = self.genotypes_df.apply(lambda x:
                                                                     self.in_interval(round(x['size_2']),
                                                                                      self.population_intervals[
                                                                                          x['micro_sat']]),
                                                                     axis=1)

    def calc_internal_relatedness(self):
        """Calculates internal relatedness and homozygosity by loci for each sample"""
        self.bin_data()
        self.ir_df = []
        for sample_name in self.genotypes_df.sample_name.unique():
            l_df = self.genotypes_df[self.genotypes_df.sample_name == sample_name]
            e_homo = 0
            e_hetero = 0
            num_hom_loci = 0
            num_loci = 0
            l_df.fillna(0.0, inplace=True)
            for row in l_df.itertuples():
                if row.binned_size_1 != 0.0 and row.binned_size_2 != 0.0:
                    num_loci += 1
                    if row.binned_size_1 == row.binned_size_2:
                        num_hom_loci += 1
                        e_homo += self.population_size_frequencies[row.micro_sat][row.binned_size_1]
                        e_homo += self.population_size_frequencies[row.micro_sat][row.binned_size_2]
                    else:
                        e_hetero += self.population_size_frequencies[row.micro_sat][row.binned_size_1]
                        e_hetero += self.population_size_frequencies[row.micro_sat][row.binned_size_2]

            score = e_homo / (e_hetero + e_homo)
            ir = (2 * num_hom_loci - (e_homo + e_hetero)) / (2 * num_loci - (e_homo + e_hetero))
            self.ir_df.append({'sample_name': sample_name, 'homozygosity_by_loci': score,
                               'internal_relatedness': ir, 'num_hom_loci': num_hom_loci,
                               'total_loci': num_loci})
        self.ir_df = pd.DataFrame(self.ir_df)
        self.ir_df = self.ir_df[['sample_name', 'homozygosity_by_loci',
                                 'internal_relatedness', 'num_hom_loci', 'total_loci']]
        self.ir_df = self.ir_df.sort_values('internal_relatedness', ascending=False)
        self.ir_df.reset_index(inplace=True, drop=True)


if __name__ == '__main__':
    print('')
