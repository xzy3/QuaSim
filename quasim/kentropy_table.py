#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 3/17/2017
#
######################

import argparse
import collections
import sys

import numpy as np
from Bio import (SeqIO)

from entropy import Profile

DEFAULT_K_RANGE=range(1,15)
DEFAULT_N=500

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
    return np.array([mean(ents) for ents in entropies]) # [variance(ents) for ents in entropies]

def join_tab_array(array):
    return "\t".join("%.3f" % a for a in array)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-k", dest="k", nargs='+', type=int, default=DEFAULT_K_RANGE)
    parser.add_argument("-n", dest="n", type=int, default=DEFAULT_N)
    parser.add_argument("-p", dest="pattern", type=str, default="")
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    fasta = list(SeqIO.parse(args.input, 'fasta'))

    ks = args.k
    patt = args.pattern
    n = args.n
    fasta = fasta[:n] if n <= len(fasta) else fasta
    actual =  get_k_entropies(fasta, ks)

    p = Profile()
    p.load_from_fasta(fasta)
    p.randomize_seqs()
    random_fasta = [s.seq for s in p.seqs]
    randomized = get_k_entropies(random_fasta, ks)

    print "\t%s" % join_tab_array(actual/randomized)
