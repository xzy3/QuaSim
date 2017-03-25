#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 3/17/2017
#
######################

import argparse
import sys
from operator import itemgetter

import collections
import numpy as np
import random
from Bio.Seq import Seq
from Bio import (AlignIO, SeqRecord, Seq, SeqIO)

from entropy import Profile

DEFAULT_K_RANGE=range(1,15)

def mean(values):
    assert isinstance(values, collections.Iterable)
    return np.mean(values)

def variance(values):
    assert isinstance(values, collections.Iterable)
    return np.var(values)


def get_k_entropies(fasta, ks):
    profiles = [Profile() for _ in ks]
    for k in ks:
        profiles[k - 1].load_from_fasta(fasta, k, patt)
    entropies = [profiles[k - 1].get_positional_entropy() for k in ks]
    print "\t".join("%.3f(%.3f)" % (mean(ents), variance(ents)) for ents in entropies)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-k", dest="k", nargs='+', type=int, default=DEFAULT_K_RANGE)
    parser.add_argument("-p", dest="pattern", type=str, default="")
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    fasta = list(SeqIO.parse(args.input, 'fasta'))

    ks = args.k
    patt = args.pattern
    get_k_entropies(fasta, ks)

    def get_shuffled_sequence(fasta, i):
        l = len(fasta[i].seq)
        seqs_array = range(len(fasta))
        str_seq = "".join(fasta[random.choice(seqs_array)].seq[i] for i in range(l))
        return SeqRecord.SeqRecord(Seq.Seq(str_seq, alphabet=Seq.Alphabet.SingleLetterAlphabet()), id=fasta[i].id)
    random_fasta = [get_shuffled_sequence(fasta, i) for i in range(len(fasta))]
    get_k_entropies(random_fasta, ks)