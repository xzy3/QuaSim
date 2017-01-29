from Bio import (SeqIO, SeqRecord, Seq)

import argparse
import sys
import math
import re
from collections import defaultdict


class Profile:

    def __init__(self, k = 1, l = 0):
        self.k = k
        self.L = l
        self.__profile__(k, l)

    def __profile__(self, k, l):
        if l > k:
            self.profile = [defaultdict(float) for _ in xrange(l - k + 1)]

    def feed_sequence(self, seq, freq = 1.0):
        assert isinstance(seq, str)
        assert isinstance(freq, float)
        seq = seq.upper()
        rng = xrange(self.L - self.k + 1)
        for i in rng:
            ss = seq[i:i+self.k]
            self.profile[i][ss] += freq

    def normalize(self):
        for e in self.profile:
            s = sum(e.values())
            for n in e:
                e[n] /= s

    def get_positional_entropy(self):
        return [sum(map(Profile.entropy_term, e.values())) for e in self.profile]

    def load_from_fasta(self, fasta, k = 1, patt = ""):
        L = -1

        for f in fasta:
            l = len(f.seq)
            if L == -1 and l > 0:
                L = l
                self.k = k
                self.L = L
                self.__profile__(k, L)
            if L != l:
                sys.stderr.write("Lengths are inconsistent %i and %i\n" % (L, l))
                sys.exit(-1)
            freq = get_count(f, patt)
            self.feed_sequence(str(f.seq), freq)

        # Normalize profile
        self.normalize()

    @staticmethod
    def entropy_term(p):
        return 0 if p == 0 else -p * math.log(p)

def get_count(fas, patt):
    if patt == "": return 1.0
    m = re.search(patt, fas.description)
    if m is not None:
        return float(m.group(0))
    return 1.0

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-k", dest="k", type=int, default=1)
    parser.add_argument("-p", dest="pattern", type=str, default="")
    args = parser.parse_args()

    fasta = SeqIO.parse(args.input, 'fasta')

    k = args.k
    patt = args.pattern

    profile = Profile()
    profile.load_from_fasta(fasta, k, patt)

    print '\n'.join(map(str, profile.get_positional_entropy()))