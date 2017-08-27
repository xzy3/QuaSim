#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 3/15/2017
#
######################
import argparse
import sys

from Bio import SeqIO
from quasim import Profile


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True,
                        help="Input fasta file. It must have all sequences of the same"
                             " length (aligned fasta) and it is required!")
    parser.add_argument("-p", dest="pattern", type=str, default="",
                        help="RegExp pattern to parse count from fasta file. Typically"
                             " it is [0-9]*$ which is number at the end of id, however"
                             " the default is \"\" which means frequency is always 1")
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout,
                        help="Output of sequence randomizer. It produces the same amount "
                             "of sequences as it is given in input file with alleles shuffled "
                             "randomly for each position. Default: stdout"
                             )
    args = parser.parse_args()

    fasta = SeqIO.parse(args.input, 'fasta')

    profile = Profile()
    profile.load_from_fasta(fasta)

    s = profile.randomize_seqs()

    profile.save_to_fasta(args.output)
