from Bio import (SeqIO, SeqRecord, Seq)

import argparse
import sys
import math
import re
import numpy
from collections import Counter

PATTERN="[0-9]*$"
FASTA='fasta'

def get_count(fas, patt):
    if patt == "": return 1.0
    m = re.search(patt, fas.description)
    if m is not None:
        return float(m.group(0))
    return 1.0

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-o", dest="output", type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    fasta = list(SeqIO.parse(args.input, FASTA))

    if len(fasta) <= 2: sys.exit(0)

    lengths = Counter(len(x) for x in fasta)

    l = lengths.most_common(1)[0][0]

    sys.stderr.write( "Initial count: %i\t" % len(fasta))
    fasta = [x for x in fasta if len(x) == l and get_count(x, PATTERN) > 2]
    sys.stderr.write("After filtering: %i\n" % len(fasta))

    SeqIO.write(fasta, args.output, FASTA)