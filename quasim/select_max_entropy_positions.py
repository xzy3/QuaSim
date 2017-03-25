#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 12/15/2016
#
######################

import argparse
import sys
from operator import itemgetter

from Bio import SeqIO
from Bio.Seq import Seq

from entropy import Profile
from quasim import PATTERN

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument("-n", dest="n", type=int, default=60)
    parser.add_argument("-p", dest="pattern", type=str, default=PATTERN)
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    try:
        fasta = list(SeqIO.parse(args.input, 'fasta'))
    except Exception:
        sys.exit(-1)

    if len(fasta) == 0:
        sys.exit(-1)
    n = args.n
    patt = args.pattern
    profile = Profile()
    profile.load_from_fasta(fasta, 1, patt)

    pos = sorted(map(itemgetter(0), sorted(enumerate(profile.get_positional_entropy()), key=itemgetter(1), reverse=True)[:n]))

    for f in fasta:
        f.seq = Seq("".join(map(lambda i: f.seq[i], pos)))

    SeqIO.write(fasta, args.output, "fasta")
