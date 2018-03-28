#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 5/1/2017
#
######################

import argparse
import itertools
import sys
import networkx as nx
import random

from Bio import (SeqIO)

from quasim import *


class Window:

    def __init__(self, start=0):
        self.onehopnet = nx.Graph()
        self.start = start
        self.end = -1
        self.fasta = None

    def find_window(self, fasta):
        cond = True
        L = len(fasta[0])
        self.end = self.start + 1
        while cond:
            self.onehopnet.clear()
            self.end = self.end + 1
            hapls = set(str(f.seq[self.start:self.end]) for f in fasta)
            self.onehopnet.add_nodes_from(hapls)
            self.onehopnet.add_edges_from(hh for hh in itertools.combinations(hapls, 2) if dist(hh[0], hh[1]) <= 1)
            # for hh in itertools.combinations(hapls, 2):
            #     if Window.dist(hh[0], hh[1]) <= 1:
            #         self.onehopnet.add_edge(hh[0], hh[1])
            cond = nx.is_connected(self.onehopnet) and self.end < L
        if self.end > L: self.end = L

    def find_window2(self, fasta):
        self.fasta = fasta
        L = len(fasta[0])
        self.end = self.start + 1


        self.end = self.end + 1
        self.__rebuild_onehopnet()

        if not nx.is_connected(self.onehopnet) or not self.end < L:
            self.end = self.start + 1
            self.__rebuild_onehopnet()
            sys.stderr.write("(%i-%i) pair of unlinked positions\n" % (self.start, self.end))
        if self.end > L:
            self.end = L
            self.__rebuild_onehopnet()

    def __rebuild_onehopnet(self):
        self.onehopnet.clear()
        hapls = set(str(f.seq[self.start:self.end]) for f in self.fasta)
        self.onehopnet.add_nodes_from(hapls)
        self.onehopnet.add_edges_from(hh for hh in itertools.combinations(hapls, 2) if dist(hh[0], hh[1]) <= 1)


    def __str__(self):
        return "%i-%i" % (self.start, self.end)

def dist(s1, s2):
    L = len(s1)
    assert len(s1) == len(s2)
    return sum(1 for i in xrange(L) if s1[i] != s2[i])

class Epistasis:

    def __init__(self, fasta = None):
        self.windows = []
        if fasta is not None:
            self.build_epistasis_windows(fasta)

    def build_epistasis_windows(self, fasta):
        cur_end = 0
        while cur_end < len(fasta[0].seq):
            w = Window()
            w.start = cur_end
            w.find_window2(fasta)
            cur_end = w.end
            self.windows.append(w)

    def from_sequence(self, str_seq):
        return [str_seq[w.start:w.end] for w in self.windows]

    def get_random_sequence(self):
        return ''.join(random.choice(w.onehopnet.nodes()) for w in self.windows)

class EpiSeq:

    def __init__(self, epistasis, str_seq):
        self.epistasis = epistasis
        self.seq = epistasis.from_sequence(str_seq)

    def __str__(self):
        return ''.join(self.seq)

    def mutate_position(self, pos):
        i , window = next(((i, w) for (i, w) in enumerate(self.epistasis.windows) if w.start <= pos < w.end), None)
        if window is None: return
        self.seq[i] = random.choice(window.onehopnet.neighbors(self.seq[i]))

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()


    fasta = list(SeqIO.parse(args.input, 'fasta'))
    epi = Epistasis(fasta)

    seq = EpiSeq(epi, epi.get_random_sequence())

    print seq
    for i in range(2):
        seq.mutate_position(i)
    print seq
    # for i in range(10):
    #     print get_sequence(ws)
    print "%s\t%i\t%s" % (args.input.name, len(epi.windows), '\t'.join(str(w.end-w.start) for w in epi.windows))