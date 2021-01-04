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
from collections import Counter
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


def read_secquence(amplicon_file):
    """Generator reading sequences contained in the file.
      :Parameters:
         fastq_file : Path to the file
    """
    with gzip.open(amplicon_file, 'rt') as f:
        line = f.readline()
        while True:
            if len(line) ==0 :
                break

            if line[0] == ">":
                seq = ""
                line2 = f.readline()
                while line2[0] != ">":
                    seq += line2[0:-1]
                    line2 = f.readline()
                    if len(line2) == 0:
                        break
                line = line2
                yield seq

def read_fasta(amplicon_file, minseqlen):
    sequencer = read_secquence(amplicon_file)
    for seq in sequencer:
        if len(seq) >= minseqlen:
            yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    sequences = read_fasta(amplicon_file, minseqlen)
    dic = {}
    for i,seq in enumerate(sequences):
        try:
            dic[seq]+=1
        except Exception:
            dic[seq]=1

    keys  = sorted(dic, key=dic.get, reverse=True)
    for key in keys:
        if dic[key] >= mincount:
            yield [key,dic[key]]

def get_chunks(sequence, chunk_size):
    l = [None,None,None,None]
    if len(sequence) >= chunk_size*4:
        for i in range(len(l)):
            l[i] = sequence[0+(i*chunk_size):chunk_size+(i*chunk_size)]
        return l

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(read, kmer_size):
    """Generator cuting kmer contained in the sequence.
      :Parameters:
         read : sequence
         kmer_size : size of the kmer
    """
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]

def get_identity(alignment_list):
    #al = nw.global_align(alignment_list[0], alignment_list[1], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
    #print(al)
    countDiff = 0
    tot = 0
    for al in alignment_list:
        countDiff+=al.count("-")
        tot = len(al)
    ident = (tot-(countDiff))
    return ident/tot
def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    sequences = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    for sequence in sequences:
        segments = get_chunks(sequence, chunk_size)

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


if __name__ == '__main__':
    main()
