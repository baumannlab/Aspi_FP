#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import pandas as pd
import numpy as np
import os
from itertools import combinations
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch


def ret_dict_attributes(x):
    '''
    >>> test_attribute ='ID=ATIG_00015626;Name=ATIG_00015626;\
                         Alias=maker-Scpiz6a_147-exonerate_protein2genome-gene-156.4;\
                         Note=Similar to ZHX1: Zinc fingers and homeoboxes protein 1 (Pongo pygmaeus);\
                         Dbxref=Gene3D:G3DSA:1.10.10.60...;\
                         Ontology_term=GO:0003677,GO:0046872;'

    >>> ret_dict_attributes(test_attribute)
    {'Alias': 'maker-Scpiz6a_147-exonerate_protein2genome-gene-156.4',
     'Dbxref': 'Gene3D:G3DSA:1.10.10.60...
     'ID': 'ATIG_00015626',
     'Name': 'ATIG_00015626',
     'Note': 'Similar to ZHX1: Zinc fingers and homeoboxes protein 1 (Pongo pygmaeus)',
     'Ontology_term': 'GO:0003677,GO:0046872'}
    '''
    y = [i for i in str(x).split(';') if i]
    z = dict([i.split('=') for i in y])

    return z


def ret_attribute(x, attribute):
    '''
    initializes a attribute dictionary using ret_dict_attributes from a gff3 attribute column.
    returns a specific dictionary entry.

    >>> test_attribute ='ID=ATIG_00015626;Name=ATIG_00015626;\
                         Alias=maker-Scpiz6a_147-exonerate_protein2genome-gene-156.4;\
                         Note=Similar to ZHX1: Zinc fingers and homeoboxes protein 1 (Pongo pygmaeus);\
                         Dbxref=Gene3D:G3DSA:1.10.10.60...;\
                         Ontology_term=GO:0003677,GO:0046872;'

    >>> ret_attribute(x,'Note')
    'Similar to ZHX1: Zinc fingers and homeoboxes protein 1 (Pongo pygmaeus)'

    '''
    dict_attributes = ret_dict_attributes(x)
    if attribute in dict_attributes:
        return dict_attributes[attribute]
    else:
        return None


def ret_tigris_gene_symbol(x):
    '''
    saves 'Note' value from attribute dictionary from ret_attribute.
    returns the gene symbol within the 'Note' value.
    >>> test_attribute ='ID=ATIG_00015626;Name=ATIG_00015626;\
                         Alias=maker-Scpiz6a_147-exonerate_protein2genome-gene-156.4;\
                         Note=Similar to ZHX1: Zinc fingers and homeoboxes protein 1 (Pongo pygmaeus);\
                         Dbxref=Gene3D:G3DSA:1.10.10.60...;\
                         Ontology_term=GO:0003677,GO:0046872;'

    >>> ret_gene_symbol(x)
    'ZHX1'
    '''
    note = ret_attribute(x, 'Note')
    if note:
        try:
            symbol = note.split('Similar to ')[1].split(':')[0]
        except IndexError:
            return note
        return symbol
    else:
        return None


def ret_tigris_gene_id(x):
    try:
        eid = ret_attribute(x, 'ID')
    except ValueError:
        eid = 'none'
    return eid


def ret_ensembl_gene_symbol(x):
    try:
        symbol = ret_attribute(x, 'Name')
    except ValueError:
        try:
            symbol = ret_attribute(x, 'description')
        except ValueError:
            symbol = ret_attribute(x, 'gene_id')

    return symbol


def ret_ensembl_gene_id(x):
    try:
        eid = ret_attribute(x, 'ID')
    except ValueError:
        eid = 'none'
    return eid

def ret_tigris_gene_name(x):
    try:
        eid = ret_attribute(x, 'Name')
    except ValueError:
        eid = 'none'
    return eid



def gff3_to_attribute_df(gff_path, ensembl=True):
    '''
    This function takes the path to a gff file as input and stores only the
    gene features wich contain the string hox or homeobox in a pandas df.
    It also creats a new column containing only the gene symbol. Defaults are set up
    for ensemble GFF3

    '''
    gff_df = pd.read_csv(gff_path, sep='\t', header=None, names=[
        'seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    gff_df = gff_df.dropna()
    gff_df.reset_index()
    gene_attributes = []
    gene_attribute = None
    if ensembl:
        feature_types = ['gene', 'chromosome', 'supercontig', "miRNA_gene",
                         "mt_gene", "pseudogene", "rRNA_gene", "snoRNA_gene",
                         "snRNA_gene"]
        for row in gff_df.itertuples():
            if row[3] in feature_types:
                gene_attribute = str(row[-1])
            gene_attributes.append(gene_attribute)
    else:
        for row in gff_df.itertuples():
            if row[3] == 'gene':
                gene_attribute = str(row[-1])

            gene_attributes.append(gene_attribute)

    print(len(gene_attributes), len(gff_df))
    gff_df['gene_attributes'] = pd.Series(gene_attributes, index=gff_df.index)

    return gff_df


def string_subset_gff3_df(gff_df, subset_regex):
    gff_df = gff_df[gff_df['gene_attributes'].str.contains(
        subset_regex, regex=True, case=False)]

    return gff_df


def add_gene_symbol(gff_df, ensembl=False):
    if ensembl:
        gff_df['gene_symbol'] = gff_df.apply(
            lambda x: ret_ensembl_gene_symbol(x['gene_attributes']), axis=1)

    else:
        gff_df['gene_symbol'] = gff_df.apply(
            lambda x: ret_tigris_gene_symbol(x['gene_attributes']), axis=1)
    return gff_df


def add_gene_id(gff_df, ensembl=False):
    if ensembl:
        gff_df['gene_id'] = gff_df.apply(
            lambda x: ret_ensembl_gene_id(x['gene_attributes']), axis=1)

    else:
        gff_df['gene_id'] = gff_df.apply(
            lambda x: ret_tigris_gene_id(x['gene_attributes']), axis=1)

    return gff_df


def add_gene_name(gff_df, ensembl=False):
    if ensembl:
        gff_df['name'] = gff_df.apply(
            lambda x: ret_ensembl_gene_symbol(x['gene_attributes']), axis=1)

    else:
        gff_df['name'] = gff_df.apply(
            lambda x: ret_tigris_gene_name(x['gene_attributes']), axis=1)

    return gff_df


def ret_gene_df(gff_path, subset_regex, ensembl=False):
    gff_df = gff3_to_attribute_df(gff_path, ensembl=ensembl)
    print(len(gff_df))
    gff_df = string_subset_gff3_df(gff_df, subset_regex)
    print(len(gff_df))
    if not gff_df.empty:
        gff_df = add_gene_symbol(gff_df, ensembl=ensembl)
        gff_df = add_gene_id(gff_df, ensembl=ensembl)
        gff_df = add_gene_name(gff_df, ensembl=ensembl)
    else:
        print('empty GFF3!')
    col_order = ['seqid', 'source', 'feature', 'start', 'end',
                 'score', 'strand', 'phase', 'gene_symbol', 'name','gene_id',
                 'attributes', 'gene_attributes']
    gff_df = gff_df[col_order]
    
    return gff_df


def check_gene_overlaps(coords):
    change = False
    for coord1, coord2 in combinations(coords, 2):
        increase = False
        a = coord1[1][0] - 1000
        b = coord1[1][1] + 1000
        c = coord2[1][0] - 1000
        d = coord2[1][1] + 1000
        if a < c and b > d:
            increase = True
        if a > c and a < d:
            increase = True
        if b > c and b < d:
            increase = True
        if a > c and b < d:
            increase = True
        if a < c and b > c:
            increase = True
        if a < d and b > d:
            increase = True
        if increase and coord1[1][2] == coord2[1][2]:
            change = True
            coord1[1][2] += 2

    return coords, change


def gather_exons_coords(gene_gff):
    exons = gene_gff[gene_gff.feature == 'exon']
    exon_coords = []
    for index, row in exons.iterrows():
        coord = [row['start'], row['end']]
        exon_coords.append(coord)
    start = min((x[0] for x in exon_coords))
    end = max((x[1] for x in exon_coords))
    gene_coords = [start, end, 1, row['strand'], row['gene_symbol']]

    return [exon_coords, gene_coords]


def coord_to_square(coord, level, color):
    start = coord[0]
    end = coord[1]
    polygon = PolygonPatch(
        Polygon(
            [(start, level - 0.3),
             (start, level + 0.3),
             (end, level + 0.3),
             (end, level - 0.3)]
        ), color=color
    )
    return polygon


def cluster_gff_genes(gff, window_size=1500000):
    clusters = []
    for seq_name in list(set(gff['seqid'])):
        scaffold_sub = gff[gff.seqid == seq_name]
        cluster_sub = pd.DataFrame()
        cluster_lim = scaffold_sub['end'].iloc[0] + window_size
        for index, row in scaffold_sub.iterrows():
            if row['start'] < cluster_lim:
                cluster_sub = cluster_sub.append(row)
            else:
                clusters.append(cluster_sub)
                cluster_sub = pd.DataFrame().append(row)
                cluster_lim = row['end'] + window_size

        clusters.append(cluster_sub)

    return clusters


def ret_largest_gene_cluster(gff3_path, subset_regex, ensembl=True, species=''):
    gff_gene_df = ret_gene_df(gff3_path, subset_regex, ensembl=ensembl)

    if not gff_gene_df.empty:
        clusters = cluster_gff_genes(gff_gene_df)
        clusters = sorted(clusters, key=lambda x: len(x), reverse=True)
        largest_cluster = clusters[0]

        return largest_cluster
    else:
        return None


def ret_sorted_cluster(gff3_path, subset_regex, ensembl=True, species=''):
    gff_gene_df = ret_gene_df(gff3_path, subset_regex, ensembl=ensembl)
    if not gff_gene_df.empty:
        clusters = cluster_gff_genes(gff_gene_df)
        clusters = sorted(clusters, key=lambda x: len(x), reverse=True)
    else:
        clusters = None
    return clusters


def plot_exons(gff_df, species='', chromosome='', subset_regex=''):
    change = True
    gene_coords = []

    for gene in gff_df.gene_id.unique():
        gene_gff = gff_df[gff_df.gene_id == gene]
        gene_coord = gather_exons_coords(gene_gff)
        gene_coords.append(gene_coord)
    gene_coords = sorted(gene_coords,
                         key=lambda x: x[1][1] - x[1][0],
                         reverse=False)
    while change == True:
        gene_coords, change = check_gene_overlaps(gene_coords)

    plt.style.use('seaborn-whitegrid')
    fig, ax = plt.subplots(ncols = 1, nrows=1, figsize=(12,2.25))
    window_start = min((x[1][0] for x in gene_coords))
    window_end = max((x[1][1] for x in gene_coords))
    window_buffer = (window_end - window_start) // 50
    ax.set_xlim(window_start - window_buffer, window_end + window_buffer)
    ax.set_ylim(0.0, 2)
    levels = []
    for gene_coord in gene_coords:
        exon_coords = gene_coord[0]
        gene_coord = gene_coord[1]
        level = gene_coord[2]
        strand = gene_coord[3]

        if strand == '-':
            symbol = '<-'
            color = 'firebrick'
        else:
            color = 'midnightblue'
            symbol = '->'

        ax.annotate(gene_coord[4], xy=(gene_coord[0] + ((gene_coord[1] - gene_coord[0]) / 3.0),
                                       level + 0.6), ha='center', fontsize=10)
        levels.append(level)
        squares = [coord_to_square(coord=coord, level=level, color=color)
                   for coord in exon_coords]
        for square in squares:
            ax.add_patch(square)

        intron_coord = range(int(gene_coord[0]), int(gene_coord[1]), 1000)
        intron_coord.append(gene_coord[1])
        intron_coord[0] = exon_coords[0][1]  
        intron_coord[-1] = exon_coords[-1][0]  
        ax.plot(intron_coord, [level] * len(intron_coord),
                symbol, zorder=0, color=color, lw=0.5, markersize=4, markeredgewidth=0.0)

    ax.set_ylim(0.0, max(levels) + 1)
    ax.yaxis.set_visible(False)
    x_ticks = np.arange(window_start, window_end,
                        ((window_end - window_start) // 15))
    x_ticks_labels = [str(round(float(x) / 1000000.0, 2)) for x in x_ticks]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks_labels,rotation=45)
    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.set_xlabel('Position(Mb)', fontsize=16)
    fig.suptitle('%s %s %s: %s..%s' %
                 (subset_regex.upper(), species, chromosome, window_start, window_end), fontsize=18)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    #fig.savefig('../fig/%s_%s_%s_%s_%s.png' %
    #            (subset_regex, species, chromosome, window_start, window_end))

    return fig, ax

def plot_cluster_series(gff3_paths, subset_regex):
    
    for gff3_path in gff3_paths:
        species = os.path.basename(gff3_path).replace('.gff3','')
        if 'tigris' in gff3_path:
            clusters = ret_sorted_cluster(gff3_path, subset_regex, ensembl=False)
        else:
            clusters = ret_sorted_cluster(gff3_path, subset_regex, ensembl=True)
        if len(clusters) > 0:
            for cluster in clusters:
                species = os.path.basename(gff3_path).replace('.gff3','')
                chromosome = cluster['seqid'].iloc[0]
                plot_exons(cluster, species=species, chromosome = chromosome, subset_regex=subset_regex)
        else:
            print('No genes match regex = %s, for %s'% (subset_regex, species))

def plot_multi_cluster_series(gff3_paths, subset_regexs):
    for subset_regex in subset_regexs:
        plot_cluster_series(gff3_paths, subset_regex)


def plot_single_cluster_series(gff3_paths, subset_regex, species, path_targ='tigris'):
    for gff3_path in gff3_paths:
        #species = os.path.basename(gff3_path).replace('.gff3','')
        if path_targ in gff3_path:
            cluster = ret_largest_gene_cluster(gff3_path, subset_regex, ensembl=False)
        else:
            cluster = ret_largest_gene_cluster(gff3_path, subset_regex, ensembl=True)
        try:
            if not cluster.empty:
                chromosome = cluster['seqid'].iloc[0]
                return plot_exons(cluster, species=species, chromosome = chromosome, subset_regex=subset_regex)
            else:
                print('No genes match regex = %s, for %s'% (subset_regex, species))
        except AttributeError:
            print('No genes match regex = %s, for %s'% (subset_regex, species))

def plot_multi_single_cluster_series(gff3_paths, subset_regexs, **kwargs):
    for subset_regex in subset_regexs:
        plot_single_cluster_series(gff3_paths, subset_regex)


if __name__ == '__main__':
    test = '/home/dut/projects/human/annotation/Homo_sapiens.GRCh38.84.chr.gff3'
    tigris = '/home/dut/projects/tigris/genome_annotation/annotations_version1/lizard_23Jun2015_piz6a_uppercase.maker.output/functional_annotation3/a_tigris_1.renamed.putative_function.iprscan.gff'
    gff_df = ret_gene_df(test, 'hoxa', ensembl=True)
    hoxa_clusters = cluster_gff_genes(gff_df)
    hoxa_clusters = sorted(hoxa_clusters, key=lambda x: len(x), reverse=True)
    test_cluster = hoxa_clusters[0]
    # print(test_cluster)

    plot_exons(test_cluster)
