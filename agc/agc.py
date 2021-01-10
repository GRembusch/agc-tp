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
    sequence = ''
    for line in myfile:
        cut_line = str(line)[2:-3] #remove b'...\n'
        if cut_line == '':
            if len(sequence) >= minseqlen:
                yield sequence
            break
        if cut_line[0] == ">":
            if len(sequence) >= minseqlen:
                yield sequence
            sequence = ''
        else:
            sequence += cut_line


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
    for key in sorted_dict:
        if seq_dict[key] >= mincount:
            yield key, seq_dict[key]

def get_chunks(sequence, chunk_size):
    """Takes a sequence and returns chunks of given size.
      :Parameters:
          sequence: The sequence to cut.
          chunk_size: Size of a chunk.
    """
    chunks = []
    for i in range(0,len(sequence)-chunk_size,chunk_size):
        chunks.append(sequence[i:i+chunk_size])
    if len(chunks) >= 4:
        return chunks
    else:
        return None

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """Generator that extracts kmers from a sequence.
      :Parameters:
          sequence: The sequence to cut.
          kmer_size: Size of a kmer
    """
    for i in range(0, len(sequence) - kmer_size +1): # remove the "\n" at the end
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Add the kmers from the given sequence to the given kmer dictionary, with the sequence id.
      :Parameters:
          kmer_dict: Dictionary of kmers.
          sequence: The sequence we add to the dict.
          id_seq : Id of the sequence.
          kmer_size: The size of the kmers.
    """
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    """Get the 8 sequences from the given kmer dict that are most similar to the given sequence.
      :Parameters:
          kmer_dict: Dictionary of kmers.
          sequence: The sequence we compare to the dict.
          kmer_size: The size of the kmers.
    """
    counter = Counter([ids for kmer in cut_kmer(sequence, kmer_size)
        if kmer in kmer_dict for ids in kmer_dict[kmer]])
    mates = [mate[0] for mate in counter.most_common(8)]
    return mates

def get_identity(alignment_list):
    """Get the identity score of 2 sequences.
      :Parameters:
          alignment_list: List containing the 2 sequences to compare.
    """
    seq1 = alignment_list[0]
    seq2 = alignment_list[1]
    nb_identical = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            nb_identical += 1
    return (nb_identical / len(seq1) * 100)

def detect_chimera(perc_identity_matrix):
    """Determines if a sequence is a chimera or not.
      :Parameters:
          perc_identity_matrix: Matrix containing the identity rate between sequences.
    """
    stddev = 0
    seq1_sim = []
    seq2_sim = []
    for perc_identity in perc_identity_matrix:
        stddev += statistics.stdev(perc_identity)
        perc1 = perc_identity[0]
        perc2 = perc_identity[1]
        if perc1 not in seq1_sim:
            seq1_sim.append(perc1)
        if perc2 not in seq2_sim:
            seq2_sim.append(perc2)
    if len(seq1_sim) >= 2 or len(seq2_sim) >= 2:
        average_stddev = stddev / len(perc_identity_matrix)
        if average_stddev > 5:
            return True
    return False

def get_identity_matrix(chunks, parents, sequence_bank, chunk_size):
    """Get the identity matrix between a sequence and 2 parents.
      :Parameters:
          chunks: Chunks from the candidate sequence.
          parents: Parent sequences from the candidate sequence.
          sequence_bank: List of sequences that are not chimeras.
          chunk_size:  Size of the chunks.
    """
    perc_identity_matrix = [[] for chunk_index in range(len(chunks))]
    for parent in parents:
        parent_chunks = get_chunks(sequence_bank[parent], chunk_size)
        for index, chunk in enumerate(chunks):
            alignment = nw.global_align(chunk, parent_chunks[index])
            identity = get_identity(alignment)
            perc_identity_matrix[index].append(identity)
    return perc_identity_matrix

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Generator of non chimeric sequences.
      :Parameters:
          amplicon_file: fasta.gz file to read.
          minseqlen: Minimal sequence length (integer).
          mincount: minimum number of sequences (integer).
          chunk_size: Size of chunks.
          kmer_size: size of kmers.
    """
    sequence_generator = dereplication_fulllength(amplicon_file,minseqlen,mincount)
    kmer_dict = {}
    sequence_bank = []
    id_sequence = 0
    perc_identity_matrix = []

    for sequence, count in sequence_generator:
        chunks = get_chunks(sequence, chunk_size)[:4]
        mates = []
        parents = []
        for chunk in chunks:
            mate = search_mates(kmer_dict,chunk,kmer_size)
            if mate not in mates:
                mates.append(mate)
                parents = common(parents,mate)
        if len(parents) >= 2:
            perc_identity_matrix = get_identity_matrix(chunks, parents[:2],
                chunk_size, sequence_bank)
        if not detect_chimera(perc_identity_matrix):# try ! instead of not
            kmer_dict = get_unique_kmer(kmer_dict,sequence,id_sequence, kmer_size)
            sequence_bank.append(sequence)
            id_sequence +=1
            yield [sequence, count]

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Regroup sequences into a greedy clustering.
      :Parameters:
          amplicon_file: fasta.gz file to read.
          minseqlen: Minimal sequence length (integer).
          mincount: minimum number of sequences (integer).
          chunk_size: Size of chunks.
          kmer_size: size of kmers.
    """
    greedy_clustering = []
    for seq_and_count in chimera_removal(amplicon_file,minseqlen,mincount,chunk_size,kmer_size):
        greedy_clustering.append(seq_and_count)
    return greedy_clustering

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """Write the OTU in a file.
      :Parameters:
          OTU_list: List of OTU to write.
          output_file: Path to the file in which to write.
    """
    index = 0
    with open(output_file, "w") as myfile:
        for sequence, count in OTU_list:
            myfile.write(">OTU_"+str(index+1)+" occurrence:"+ str(count) + "\n")
            myfile.write(fill(str(sequence)))
            myfile.write("\n")
            index += 1


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon_file = args.amplicon_file
    if isfile(amplicon_file):
        otu_cluster = abundance_greedy_clustering(amplicon_file, args.minseqlen,args.mincount,
            args.chunk_size,args.kmer_size)
        write_OTU(otu_cluster, args.output_file)

if __name__ == '__main__':
    main()
