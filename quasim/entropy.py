
from Bio import (SeqIO, SeqRecord, Seq)

import random
import argparse
import sys
import math
import re
from collections import defaultdict
import numpy as np

PATTERN=""


class Variant:

    def __init__(self, seq, patt = PATTERN):
        assert isinstance(seq, SeqRecord.SeqRecord)

        self.seq = seq
        self.freq = get_count(seq, patt)
        self.id = seq.id

class Profile:

    def __init__(self, k = 1, l = 0):
        self.k = k
        self.L = l
        self.seqs = []
        self.__profile__(k, l)
        self.mut_map = []
        self.size = 0

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

    def distance_to(self, seq):
        seq = str(seq)
        assert len(seq) == len(self.profile)
        d = 0.0
        for i,nucl in enumerate(seq):
            for n in self.profile[i]:
                if n == nucl:
                    d += (1-self.profile[i][n])**2
                else:
                    d += self.profile[i][n]**2
        return math.sqrt(d)

    def std(self):
        vals = [(self.distance_to(s[0]), s[1]) for s in self.seqs]
        size = sum(v[1] for v in vals)
        m = sum(v[0]*v[1] for v in vals)
        m /= size
        sigma = .0
        for v in vals:
            sigma += ((v[0]-m)**2) * v[1]
        sigma = math.sqrt(sigma/size)
        self.size = size
        return m, sigma

    def load_from_fasta(self, fasta, k = 1, patt = PATTERN):
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
                continue
            v = Variant(f, patt)
            self.seqs.append(v)
            self.feed_sequence(str(f.seq), v.freq)

        # Normalize profile
        self.normalize()

    def save_to_fasta(self, fasta):

        SeqIO.write([s.seq for s in self.seqs], fasta, "fasta")

    def randomize_seqs(self):
        str_seqs = [str(s.seq.seq) for s in self.seqs]
        str_seqs = Profile.__randomize_seqs(str_seqs)

        for i in range(len(self.seqs)):
            self.seqs[i].seq = SeqRecord.SeqRecord(
                                    Seq.Seq(str_seqs[i], alphabet=Seq.Alphabet.SingleLetterAlphabet()),
                                    id=self.seqs[i].id, description='')

    def build_mutations_map(self):
        if self.k == 1:
            self.mut_map = [{x: [y for y in e if y != x] for x in e} for e in self.profile]

    def mutate_nucl_pos(self, nucl, pos):
        """
        Check if mut_map has mutation for a given
        position and nucleotide and return another
        nucleotide if possible
        :param nucl: Original nucleotide
        :param pos:  Position in the region
        :return: New nucleotide or None
        """
        if 0 <= pos < len(self.mut_map) and nucl in self.mut_map[pos]:
            return random.choice(self.mut_map[pos][nucl])
        print "Pos: %i Nucl: %s" % (pos, nucl)
        return None

    @staticmethod
    def __randomize_seqs(seqs):
        if len(seqs) == 0:
            return None

        L = len(seqs[0])
        N = range(len(seqs))
        rnd_seqs = [list(s) for s in seqs]
        rnd_indices = [np.random.permutation(N) for _ in range(L)]

        for i in N:
            rnd_seqs[i] = [seqs[rnd_indices[x][i]][x] for x in range(L)]

        return [''.join(x) for x in rnd_seqs]

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
    parser.add_argument("-p", dest="pattern", type=str, default=PATTERN)
    args = parser.parse_args()

    fasta = SeqIO.parse(args.input, 'fasta')

    k = args.k
    patt = args.pattern

    profile = Profile()
    profile.load_from_fasta(fasta, k, patt)

    #print '\n'.join(map(str, profile.get_positional_entropy()))
    print '\n'.join('\t'.join(map(str, e.values())) for e in profile.profile)
    profile.build_mutations_map()
    print profile.mut_map