#!/home/dut/anaconda2/bin/python
from __future__ import print_function
import numpy as np
import re
import argparse
from collections import Counter
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from IPython.display import Markdown, display, HTML


def write_two_line_fasta(file_name, a_fasta_tuples):
    fh = open(file_name, "w")
    for i in np.arange(len(a_fasta_tuples)):
        fh.write(str(a_fasta_tuples[i][0]) + "\n" +
                 str(a_fasta_tuples[i][1]) + "\n")


def print_fasta_tuple(fasta_tuples):
    for i in range(len(fasta_tuples)):
        print(fasta_tuples[i][0])
        print(fasta_tuples[i][1])


def print_fasta_dict(fasta_dict):
    for key, val in fasta_dict.items():
        print(key)
        print(val)

# ########################################################################
# ############################# Stats Functions ##########################
# ########################################################################


def calcNX(list_of_lengths, percent):
    ''' This function takes as input an list of lengths and a percentage. It
    then calculates and returns the NX (for X = 50%, this calculates the N50).
    Given a set of sequences of varying lengths, the NX length is defined as
    the length N for which X% of all bases in the sequences are in a sequence
    of length L < N.
    Example:
    list_of_lengths = [5, 10, 15, 16, 17, 20]
    p = 0.50
    N50 = calcNX(list_of_lengths, p)
    print N50
    16
    '''
    list_of_lengths = sorted(list_of_lengths)
    x_distance = sum(list_of_lengths) * percent
    running_length = 0
    i = 0
    while running_length < x_distance:
        running_length += list_of_lengths[i]
        i += 1
    NX = list_of_lengths[i - 1]
    return NX


def nxDistribution(list_of_lengths):
    nx_values = []
    nx_percents = []
    for n in range(1, 100):
        p = float(n / 100.0)
        x = calcNX(list_of_lengths, p)
        z = int(100 - n)
        nx_percents.append(z)
        nx_values.append(x)
    nx_dist = [nx_percents, nx_values]
    return nx_dist


def listMedian(list_of_lengths):
    '''
    This function take an list of lengths as input, and returns the median
    '''
    list_of_lengths = sorted(list_of_lengths)
    if len(list_of_lengths) % 2 == 0:
        a = list_of_lengths[int(len(list_of_lengths) / 2)]
        b = list_of_lengths[int((len(list_of_lengths) / 2) - 1)]
        median = (a + b) / 2
    else:
        median = list_of_lengths[int((len(list_of_lengths) - 1) / 2)]
    return median


def getKeyForMaxValue(d):
    '''
    Found this on stack overflow. Takes a dictionary and returns the key
    corresponding to the max value in the dictionary. The stack overflow
    entry contained tests showing this was the fastest method.
    http://stackoverflow.com/a/12343826
    '''
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


def lengths_greater_than(list_of_lengths, lim):
    filt_list = [length for length in list_of_lengths if length > lim]
    return filt_list


def num_lengths_greater_than(list_of_lengths, lim):
    x = len(lengths_greater_than(list_of_lengths, lim))
    return x


def perc_greater_than(list_of_lengths, lim):
    total = sum(list_of_lengths)
    greater_than = sum(lengths_greater_than(list_of_lengths, lim))
    perc = (float(greater_than) / float(total)) * 100
    return perc


# ##############################################################################
# ########################### Sequence Functions #########################
# ##############################################################################

def get_reverse_compliment(seq):
    seq = seq[::-1]
    rev_comp = []
    for char in seq:
        if char == "T":
            rev_comp.append("A")
        elif char == "A":
            rev_comp.append("T")
        elif char == "G":
            rev_comp.append("C")
        elif char == "C":
            rev_comp.append("G")
        else:
            rev_comp.append(char)

    rev_comp = "".join(rev_comp)
    return rev_comp


def translateDNA(dna):
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
    'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
    }
    protein=''
    for i in range(0,len(dna),3):
        codon=dna[i:i+3]
        if len(codon)==3:
            protein+=codontable[codon]
    return protein


def atg_translate(dna):
    atg_proteins = []
    for i in xrange(len(dna) - 3 + 1):
        codon = dna[i:i+3]
        if codon == 'ATG':
            cdna= dna[i:]
            protein = translateDNA(cdna)
            atg_proteins.append(protein)
    return atg_proteins

def atg_stop_translate(dna):
    atg_proteins = []
    for i in xrange(len(dna) - 3 + 1):
        codon = dna[i:i+3]
        if codon == 'ATG':
            protein = ''
            cdna= dna[i:]
            for i in xrange(len(cdna) - 3 + 1):
                codon = cdna[i:i+3]
                aa = translateDNA(codon)
                if aa != '-':
                    protein += aa
                else:
                    atg_proteins.append(protein)
                    break
    atg_proteins = sorted(atg_proteins,key=len)[::-1]
    return atg_proteins


def find_ortholog_translation(dna, ortho):
    forward = {i + 1:atg_stop_translate(dna[i:])[0] for i in xrange(3)}
    reverse = {-1 * i +1 : atg_stop_translate(get_reverse_compliment(dna)[i:])[0] for i in xrange(3)}
    all_frames =  forward.copy()
    all_frames.update(forward)
    alignmets = {k:pairwise2.align.localxx(ortho, v.replace('-',''),  score_only=True) for k, v in all_frames.items()}
    best_frame = getKeyForMaxValue(alignmets)

    return (best_frame, all_frames[best_frame])

def extracBlasttORF(seq, start, stop):
    if start > stop:
        seq = seq[::-1]
        start, stop = stop, start
    start = start-1
    orf = seq[start:stop]

    return orf

def baseComposition(seq):
    '''
    This funciton takes a string as input and returns the percentage of A,T,G,
    C and N that occur in the string in a dictionary
    '''
    composition = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
    total = 0
    for i in seq:
        if i in composition:
            composition[i] += 1
            total += 1
    for x in composition:
        composition[x] = composition[x] / float(total) * 100
    return composition


def ret_gc_content(seq, missing_cutoff=1.0):
    missing_cutoff = missing_cutoff * 100
    composition = baseComposition(seq)
    if composition['N'] >= missing_cutoff:
        gc = 0.0
    else:
        gc = composition['G'] + composition['C']
    return gc


def kmerize_seq(sequence, size):
    c = len(sequence)
    kmers = []
    for i in np.arange(c - size + 1):
        start = i
        stop = i + size
        kmer = sequence[start:stop]
        kmers.append(kmer)
    return kmers


def nuc_stretch(seq, Nuc):
    '''
    This function takes string as input and finds all stretches of a certain
    letter, Nuc, it then returns a list of all lengths of all stretches found
    in the string
    '''
    if Nuc in seq:
        nuc_length = [len(a) for a in re.findall(r'%c+' % (Nuc), seq)]
        return nuc_length
    else:
        return [0]


def convert_string_html(string):
    html_string = ''.join(['<span><tt>%s</tt></span>'%s for s in string])
    return html_string


def display_html_alignment(alignment,title ='', matrix = matlist.blosum62):
    html_alignment='<span ><tt><b>%s</b></tt></span><br />'% title
    all_strings =[]
    for i in range(0,len(alignment[0]),80):
        string1 = convert_string_html(alignment[0][i:i+80])
        similarity =''
        for k in range(len(alignment[0][i:i+80])):
            aa1 = alignment[0][i:i+80][k]
            aa2 = alignment[1][i:i+80][k]
            if aa1 == aa2:
                similarity += '<span style="background-color: #2ca25f" ><tt>%s</tt></span>'%"*"
            elif aa1 != "-" and aa2 != "-":
                try:
                    blosum62 = matrix[(aa1, aa2)]
                except KeyError:
                    blosum62 = matrix[(aa2, aa1)]
                if blosum62 < 0:
                    similarity += '<span style="background-color: #e5f5f9"><tt>%s</tt></span>'%"."
                else:
                    similarity += '<span style="background-color: #99d8c9"><tt>%s</tt></span>'%":"
            else:
                similarity += '<span><tt>%s</tt></span>'%"|"

        string2 = similarity
        string3 = convert_string_html(alignment[1][i:i+80])
        all_strings.append('<br />'.join([string1,string2,string3]) + '<br></br>')
    html_alignment += '<br />'.join(all_strings) +'<br />'
    display(HTML(html_alignment))


def display_pairwise_alignment(seq1, seq2, title='', matrix = matlist.blosum62, gap_open = -10, gap_extend = -0.5):
    alignment = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)[0]
    display_alignment(alignment,title,matrix)

###############################################################################
############################     Fasta Class    ###############################
###############################################################################


class Fasta_file(object):
    '''
    This is a fasta file class, it takes a fasta file as input and transforms
    it into various object types that a user may need.
    '''

    def __init__(self, path):
        self.path = path
        fasta_fh = open(self.path, "r")
        fasta_dict = {}
        for line in fasta_fh:
            line = line.strip()
            if ">" in line:
                try:
                    fasta_dict[header] = "".join(fasta_dict[header])
                except UnboundLocalError:
                    pass
                header = line.strip()
                fasta_dict[header] = []
                continue
            fasta_dict[header].append(line)
        fasta_dict[header] = "".join(fasta_dict[header])
        fasta_fh.close()

        self.stats = None
        self.bc = None
        self.fasta_dict = fasta_dict
        self.fasta_tuples = fasta_dict.items()

    def unmask(self):
        self.fasta_tuples = [(x[0], x[1].upper()) for x in self.fasta_tuples]
        self.fasta_dict = dict(self.fasta_tuples)

    def sort_fasta_tuples(self):
        self.fasta_tuples = sorted(self.fasta_tuples,
                                   key=lambda x: len(x[1]),
                                   reverse=True)

    def size_filter_tuples(self, min_length):
        min_length = int(min_length)
        size_filt_tuples = []
        for i in np.arange(len(self.fasta_tuples)):
            if len(self.fasta_tuples[i][1]) > min_length:
                size_filt_tuples.append(self.fasta_tuples[i])

        return size_filt_tuples

    def size_filter_dict(self, min_length):
        filtered = self.size_filter_tuples(min_length)
        size_filt_dict = dict(filtered)
        return size_filt_dict

    def size_filt_fasta(self, min_length):
        self.fasta_tuples = self.size_filter_tuples(min_length)
        self.fasta_dict = dict(self.fasta_tuples)

    def subset_fasta_tuples(self, id_list):
        id_list = [">" + i.replace(">", "") for i in id_list]
        subsetted_fasta = []
        for seq_id in id_list:
            try:
                subsetted_fasta.append((seq_id, self.fasta_dict[seq_id]))
            except KeyError:
                print(seq_id + " Not in fasta_dict")
        return subsetted_fasta

    def subset_fasta_dict(self, id_list):
        subsetted_fasta = dict(self.subset_fasta_tuples(id_list))
        return subsetted_fasta

    def regex_subset_fasta_dict(self, regex):
        subsetted_fasta = {key:val for key,val in self.fasta_dict.items() if re.search(regex,key,re.IGNORECASE)}
        return subsetted_fasta

    def extract_sequence(self, scaffold, start, end, name=''):
        '''
        Extract a sequence region from a sequence
        Takes a header key, start position, stop position and optionally a name
        Uses 1 start counting (first nucleotide is 1)
        '''
        extract_id = '%s:%d..%d %s' % (scaffold, start, end, name)
        extract_seq = self.fasta_dict[scaffold][start - 1:end]
        return (extract_id, extract_seq)

    def extract_multi_seqs(self,list_of_targets):
        extracted_seqs = []
        for target in list_of_targets:
            extracted_seqs.append(self.extract_sequence(*target))
        return extracted_seqs

    def preview_fasta(self):
        for i in range(5):
            print(self.fasta_tuples[i])

    def strip_headers(self):
        stripped_fasta_tuples = [(x[0].split()[0], x[1])
                                 for x in self.fasta_tuples]
        self.fasta_tuples = stripped_fasta_tuples
        self.fasta_dict = dict(stripped_fasta_tuples)
        return self

    def add_meta_header(self, meta_dict):
        self.strip_headers()
        temp = {}
        for key in self.fasta_dict:
            if key.strip(">") in meta_dict:
                temp[key + " " + str(meta_dict[key.strip(">")])
                     ] = self.fasta_dict[key]
            else:
                temp[key] = self.fasta_dict[key]

        self.fasta_dict = temp
        self.fasta_tuples = temp.items()

    def make_length_list(self):
        self.length_list = sorted([len(x[1]) for x in self.fasta_tuples])

    def make_genome_string(self):
        self.genome_string = "".join(self.fasta_dict.values())

    def get_bc(self):
        self.make_genome_string()
        self.bc = baseComposition(self.genome_string)

    def get_stats(self):
        self.make_length_list()
        self.stats = {"max_scaffold": max(self.length_list),
                      "min_scaffold": min(self.length_list),
                      "mean_scaffold": np.mean(self.length_list),
                      "median_scaffold": listMedian(self.length_list),
                      "N50": calcNX(self.length_list, 0.5),
                      "N90": calcNX(self.length_list, 0.1),
                      "scaffolds_greater_than_10kb":
                      num_lengths_greater_than(self.length_list, 10000),
                      "scaffolds_greater_than_100kb":
                      num_lengths_greater_than(self.length_list, 100000),
                      "scaffolds_greater_than_1Mb":
                      num_lengths_greater_than(self.length_list, 1000000),
                      "scaffolds_greater_than_10Mb":
                      num_lengths_greater_than(self.length_list, 10000000),
                      "scaffolds_greater_than_100Mb":
                      num_lengths_greater_than(self.length_list, 100000000),
                      "perc_genome_covered_by_scaffolds_greater_than_10kb":
                      perc_greater_than(self.length_list, 10000),
                      "perc_genome_covered_by_scaffolds_greater_than_100kb":
                      perc_greater_than(self.length_list, 100000),
                      "perc_genome_covered_by_scaffolds_greater_than_1mb":
                      perc_greater_than(self.length_list, 1000000),
                      "number of scaffolds": len(self.length_list),
                      "total bases": sum(self.length_list)
                      }
        for key, val in self.stats.items():
            print(key, val)

    def get_bc_stats(self):
        if not self.stats:
            self.get_stats()
        if not self.bc:
            self.get_bc()
        self.stats["perc_A"] = self.bc["A"]
        self.stats["perc_T"] = self.bc["T"]
        self.stats["perc_G"] = self.bc["G"]
        self.stats["perc_C"] = self.bc["C"]
        self.stats["perc_N"] = self.bc["N"]

    def ret_n_len_dist(self):
        dist = []
        for seq in self.fasta_dict.values():
            dist.extend(nuc_stretch(seq))
        dist = Counter(dist)

        return dist

    def write_n_len_dist(self):
        filename = self.path.replace(".fasta", ".n_len_dist.tsv")
        filename = filename.replace(".fa", ".n_len_dist.tsv")
        fho = open(filename)
        dist = sorted(self.ret_n_len_dist().items(),
                      key=lambda x: x[1], reverse=True)
        for x in dist:
            fho.write('%s\t%s\n' % (x[0], x[1]))
        fho.close()

    def write_stats(self):
        filename = self.path.replace(".fasta", ".statistics.tsv")
        filename = filename.replace(".fa", ".statistics.tsv")
        fh = open(filename, "w")
        if not self.stats:
            self.get_stats()
        for key, val in self.stats.items():
            print(key, val)
            fh.write("%s\t%.2f\n" % (key, val))
        fh.close()

    def return_nx_dist(self):
        self.make_length_list()
        nx_dist = nxDistribution(self.length_list)
        return nx_dist

    def get_sliding_window_gc_data(self, window_size=10000, missing_cutoff=0.20):
        self.make_genome_string()
        self.sliding_window_data_gc_data = []

        for i in np.arange(0, len(self.genome_string) - window_size, window_size):
            window = self.genome_string[i:i + window_size]
            gc = ret_gc_content(window, missing_cutoff=missing_cutoff)
            window_data = [i, i + window_size, gc]
            self.sliding_window_data_gc_data.append(window_data)

        last_window = self.genome_string[i + window_size:]
        if len(last_window) > (1 - missing_cutoff) * window_size:
            gc = ret_gc_content(last_window, missing_cutoff=missing_cutoff)
            window_data = [i + window_size, len(self.genome_string), gc]
            self.sliding_window_data_gc_data.append(window_data)

    def get_sliding_scaffold_window_gc_data(self, window_size=10000, missing_cutoff=0.20):
        self.sliding_scaffold_window_data_gc_data = []
        for header in self.fasta_dict:
            if len(self.fasta_dict[header]) < window_size + 1:
                continue

            for i in np.arange(0, len(self.fasta_dict[header]) - window_size, window_size):
                window = self.fasta_dict[header][i:i + window_size]
                gc = ret_gc_content(window, missing_cutoff=missing_cutoff)
                window_data = [header, i, i + window_size, gc]
                self.sliding_scaffold_window_data_gc_data.append(window_data)

            last_window = self.fasta_dict[header][i + window_size:]
            if len(last_window) > (1 - missing_cutoff) * window_size:
                gc = ret_gc_content(last_window, missing_cutoff=missing_cutoff)
                window_data = [header, i + window_size,
                               len(self.fasta_dict[header]), gc]
                self.sliding_scaffold_window_data_gc_data.append(window_data)

    def write_sliding_gc_data(self, window_size=10000):
        filename = self.path.replace(
            ".fasta", ".gc_%s_window.tsv" % str(window_size))
        filename = filename.replace(
            ".fa", ".gc_%s_window.tsv" % str(window_size))
        fh = open(filename, "w")
        self.get_sliding_window_gc_data(window_size=window_size)
        for data in self.sliding_window_data_gc_data:
            fh.write("%d\t%d\t%.2f\n" % (data[0], data[1], data[2]))
        fh.close()

    def write_sliding_scaffold_gc_data(self, window_size=10000):
        filename = self.path.replace(
            ".fasta", ".gc_%s_scaffold_window.tsv" % str(window_size))
        filename = filename.replace(
            ".fa", ".gc_%s_scaffold_window.tsv" % str(window_size))
        fh = open(filename, "w")
        self.get_sliding_scaffold_window_gc_data(window_size=window_size)
        for data in self.sliding_scaffold_window_data_gc_data:
            fh.write("%s\t%d\t%d\t%.2f\n" %
                     (data[0], data[1], data[2], data[3]))
        fh.close()


def run(args):
    fasta = Fasta_file(args.fasta)
    if args.unmask:
        fasta.unmask()
    if args.filt:
        fasta.size_filt_fasta(args.filt)
    if args.write:
        write_two_line_fasta(args.write, fasta.fasta_tuples)
    if args.bc_stats:
        fasta.get_bc_stats()
        fasta.write_stats()
    elif args.stats:
        fasta.write_stats()
    if args.gc_windows:
        fasta.write_sliding_gc_data()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-fasta', action='store', dest='fasta',
                        help="path to fasta_file")

    parser.add_argument('-filt', action='store', dest='filt', default=None,
                        help='only use scaffolds greater_than this ammount')

    parser.add_argument('-unmask', action='store_true', dest='unmask', default=False,
                        help='convert all lowercase letters to uppercase')

    parser.add_argument('-write_fasta', action='store', dest='write', default=False,
                        help='write an fasta file with applied filters')

    parser.add_argument('-stats', action='store_true', dest='stats',
                        help="write length stats")

    parser.add_argument('-bc_stats', action='store_true', dest='bc_stats',
                        help='write length and base composition stats')

    parser.add_argument('-gc_windows', action='store_true', dest='gc_windows',
                        help='write out 10kb GC content windows')

    args = parser.parse_args()
    print(args)
    run(args)
