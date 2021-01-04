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

__author__ = "Moreau Baptiste"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Moreau Baptiste"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Moreau Baptiste"
__email__ = "moreaubapt@eisti.eu"
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
    """
    prend deux arguments correspondant au fichier fasta.gz et à la longueur minimale
    des séquences et retourne un générateur de séquences de longueur l >= minseqlen: yield sequence.
    """
    sequencer = read_secquence(amplicon_file)
    for seq in sequencer:
        if len(seq) >= minseqlen:
            yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Prend trois arguments correspondant au fichier fasta,  la longueur minimale des séquences
    et leur comptage minimum. Elle fait appel au générateur fourni par read_fasta
    et retourne un générateur des séquences uniques ayant une occurrence O>=mincount
    ainsi que leur occurrence. Les séquences seront retournées par ordre décroissant d’occurrence:
    yield [sequence, count]
    """
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
    """
    prend une séquence et une longueur de segment l: chunk_size et retourne une liste de sous-séquences de taille l non chevauchantes. A minima 4 segments doivent être obtenus par séquence.

    """
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
    """prend un alignement (sous forme de liste) et calcule le pourcentage d’identité entre
    es deux séquences selon la formule: id = nb nucléotides identiqueslongueur de l'alignement
    """
    total = len(alignment_list[0])
    identik = 0
    for i in range(total):
        if alignment_list[0][i] == alignment_list[1][i]:
            identik +=1
    return 100*(identik/total)

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Fait appel au générateur fourni par dereplication_fulllength et
    retourne un générateur des séquences non chimérique au format:
    yield [sequence, count]
    """
    kmer_dict = {}
    perc_identity_matrix = []
    chunk_match = []
    seq_list = []
    chim_id = 0
    for i, occurence_list in enumerate(list(dereplication_fulllength(amplicon_file, minseqlen,
                                                                     mincount))):
        chim = True
        chunk_list = get_chunks(occurence_list[0], chunk_size)
        for chunk in chunk_list:
            chunk_match.append([i[0] for i in Counter([ids for kmer in cut_kmer(chunk, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)])
        com_seq = common(chunk_match[0], chunk_match[1])
        for j in range(2, len(chunk_match)):
            com_seq = common(com_seq, chunk_match[j])
        if len(com_seq) > 1:
            for k in range(len(chunk_list)):
                perc_identity_matrix.append([][k])
            for seq in com_seq[0:2]:
                seq_chunk_list = get_chunks(seq_list[seq], chunk_size)
                for l, chunk in enumerate(chunk_list):
                    perc_identity_matrix[l].append(get_identity(
                        nw.global_align(chunk, seq_chunk_list[l],
                            gap_open = -1, gap_extend = 1, matrix = "MATCH")))
            std_list = []
            flag_file = 0
            flag_similarity = 0
            for line in perc_identity_matrix:
                std_list.append(statistics.stdev([line[0], line[1]]))
                if flag_file == 0:
                    val0 = line[0]
                    val1 = line[0]
                    flag_file = 1
                else:
                    if flag_similarity == 1:
                        continue
                    if val0 != line[0] and val1 != line[1]:
                        flag_similarity = 1
                    val0 = line[0]
                    val1 = line[0]
            std_mean = statistics.mean(std_list)
            if std_mean > 5 and flag_similarity == 1:
                chim = True
            chim =  False



        else:
            chim = False
        if not chim:
            #kmer_dict = get_unique_kmer(kmer_dict, occurence_list[0], chim_id, kmer_size)
            for kmer in cut_kmer(occurence_list[0],kmer_size):
                if kmer not in kmer_dict:
                    kmer_dict[kmer] = [chim_id]
                else :
                    kmer_dict[kmer].append(chim_id)
                    kmer_dict[kmer] = list(set(kmer_dict[kmer]).union(set(kmer_dict[kmer])))
            seq_list.append(occurence_list[0])
            chim_id += 1
            yield occurence_list

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
