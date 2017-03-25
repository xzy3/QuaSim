from Bio import (SeqIO, SeqRecord, Seq)

import argparse
import sys
import math
import re
# ! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 3/15/2017
#
######################

import numpy
from collections import Counter

from quasim import FASTA, PATTERN, DEFAULT_MIN_SIZE


def get_count(fas, patt):
    if patt == "": return 1.0
    m = re.search(patt, fas.description)
    if m is not None:
        return float(m.group(0))
    return 1.0

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-n", dest='n', type=int, default=DEFAULT_MIN_SIZE)
    parser.add_argument("-o", dest="output", type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    fasta = list(SeqIO.parse(args.input, FASTA))

    if len(fasta) <= 2: sys.exit(0)

    lengths = Counter(len(x) for x in fasta)

    l = lengths.most_common(1)[0][0]

    sys.stderr.write( "Initial count: %i\t" % len(fasta))
    fasta = [x for x in fasta if len(x) == l and get_count(x, PATTERN) > 2]
    sys.stderr.write("After filtering: %i\n" % len(fasta))

    if len(fasta) < args.n:
        sys.stderr.write("Sample is too small! Skipping...\n")
    else:
        SeqIO.write(fasta, args.output, FASTA)