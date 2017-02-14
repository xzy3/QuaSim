import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq

def distance(s1, s2):
    if isinstance(s1, Seq):
        s1 = str(s1.seq)
    if isinstance(s2, Seq):
        s2 = str(s2.seq)
    assert len(s1) == len(s2)
    return sum(s != t for (s, t) in zip(s1, s2))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i1", dest='input1', type=argparse.FileType('r'), required=True)
    parser.add_argument("-i2", dest='input2', type=argparse.FileType('r'), required=True)
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    fasta1 = list(SeqIO.parse(args.input1, 'fasta'))
    fasta2 = list(SeqIO.parse(args.input2, 'fasta'))

    den = len(fasta1) * len(fasta2)

    num = sum(distance(s, t) for s in fasta1 for t in fasta2)

    args.output.write("%.3f\n" % (float(num)/den))
