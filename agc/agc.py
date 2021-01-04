#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter, defaultdict
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """Generator that reads the given file.
      :Parameters:
          amplicon_file: fasta.gz file to read
          minseqlen: Minimal sequence length (integer)
    """
    myfile = gzip.open(amplicon_file)
    for line in myfile: # >Name
        sequence = ''
        current_line = str(line)[2:-3]
        if current_line == '':
            break
        if current_line[0] == '>':
            current_line = str(next(myfile))[2:-3]
        while current_line != '' and current_line[0] != '>':
            sequence += current_line
            current_line = str(next(myfile))[2:-3]
        if len(sequence) >= minseqlen:
            yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Generator of sequences for which the number of occurences is higher than mincount
      :Parameters:
          amplicon_file: fasta.gz file to read
          minseqlen: Minimal sequence length (integer)
          mincount: minimum number of sequences (integer)
    """
    seq_dict = defaultdict(lambda: 0)
    for sequence in read_fasta(amplicon_file,minseqlen):
        seq_dict[sequence] += 1
    sorted_occurences = sorted(seq_dict.items(), key = lambda kv: kv[1], reverse=True)
    sorted_dict = dict(sorted_occurences)
    print(sorted_dict)
    for key in sorted_dict:
        if seq_dict[key] >= mincount:
            yield key, seq_dict[key]



def get_chunks(sequence, chunk_size):
    pass

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    pass

def get_identity(alignment_list):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    pass
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    for [key,value] in dereplication_fulllength(args.amplicon_file,200,3):
        print(key)


if __name__ == '__main__':
    main()
