import argparse
import sys
from operator import itemgetter

from Bio import SeqIO
from Bio.Seq import Seq

from entropy import Profile

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-n", dest="n", type=int, default=60)
    parser.add_argument("-p", dest="pattern", type=str, default="")
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    fasta = list(SeqIO.parse(args.input, 'fasta'))

    n = args.n
    patt = args.pattern
    profile = Profile()
    profile.load_from_fasta(fasta, 1, patt)

    pos = sorted(map(itemgetter(0), sorted(enumerate(profile.get_positional_entropy()), key=itemgetter(1), reverse=True)[:n]))

    for f in fasta:
        f.seq = Seq("".join(map(lambda i: f.seq[i], pos)))

    SeqIO.write(fasta, args.output, "fasta")
