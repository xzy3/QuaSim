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
DEFAULT_MAX_T=1000

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
    print "\t".join("%.3f" % (mean(ents)) for ents in entropies)
    # , variance(ents)

def find_best_t(input, t, reads, patt):
    L = len(reads)
    reads_pr = Profile()
    reads_pr.load_from_fasta(reads)
    real_mean = mean(reads_pr.get_positional_entropy())
    best_T=-1
    best_seqs=[]
    diff = 10000
    for i in range(1, t):
        sim_seqs = list(SeqIO.parse(input % i, "fasta"))
        if len(sim_seqs) > L:
            sub_seqs = random.sample(sim_seqs, L)
        else:
            sub_seqs = sim_seqs
        pr = Profile()
        pr.load_from_fasta(sub_seqs)
        sim_mean = mean(pr.get_positional_entropy())
        if abs(sim_mean - real_mean) < diff:
            diff = abs(sim_mean - real_mean)
            best_T = i
            best_seqs = sim_seqs

    return best_T, best_seqs



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=str, required=True)
    parser.add_argument("-t", dest='t', type=int, default=DEFAULT_MAX_T)
    parser.add_argument("-r", dest='reads', type=argparse.FileType('r'), required=True)
    parser.add_argument("-k", dest="k", nargs='+', type=int, default=DEFAULT_K_RANGE)
    parser.add_argument("-p", dest="pattern", type=str, default="")
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    reads = list(SeqIO.parse(args.reads, 'fasta'))
    ks = args.k
    patt = args.pattern

    T, sim_reads = find_best_t(args.input, args.t, reads, patt)

    print T

    fasta = sim_reads

    get_k_entropies(reads, ks)

    get_k_entropies(fasta, ks)

    def get_shuffled_sequence(reads, i):
        l = len(reads[i].seq)
        seqs_array = range(len(reads))
        str_seq = "".join(reads[random.choice(seqs_array)].seq[i] for i in range(l))
        return SeqRecord.SeqRecord(Seq.Seq(str_seq, alphabet=Seq.Alphabet.SingleLetterAlphabet()), id=reads[i].id)
    random_fasta = [get_shuffled_sequence(reads, i) for i in range(len(reads))]
    get_k_entropies(random_fasta, ks)