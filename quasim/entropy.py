
from Bio import (SeqIO, SeqRecord, Seq)

import random
import argparse
import sys
import math
import re
from collections import defaultdict


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
                continue
            freq = get_count(f, patt)
            self.seqs.append((str(f.seq), freq))
            self.feed_sequence(str(f.seq), freq)

        # Normalize profile
        self.normalize()

    def save_to_fasta(self, fasta):

        s = [(seq_freq[0], i) for i, seq_freq in enumerate(self.seqs)]
        seqs = [SeqRecord.SeqRecord(Seq.Seq(str_seq, alphabet=Seq.Alphabet.SingleLetterAlphabet()),
                                    id="%i_%.0f" % (i, 1.0), description='')
                    for str_seq, i in s]
        SeqIO.write(seqs, fasta, "fasta")

    def build_mutations_map(self):
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
        return None

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

    #print '\n'.join(map(str, profile.get_positional_entropy()))
    print '\n'.join('\t'.join(map(str, e.values())) for e in profile.profile)
    profile.build_mutations_map()
    # print profile.mut_map