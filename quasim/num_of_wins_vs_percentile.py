#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 7/18/2017
#
######################

import argparse
import itertools
import sys


from Bio import (SeqIO)

from quasim import *
from epistasis import *
from entropy import get_count

def get_number_of_variants_percentile(fasta, percent=0.5):
    if percent < 0 or percent > 1:
        return None
    patt = "[0-9]*$"

    counts = sorted([get_count(fas, patt) for fas in fasta], reverse=True)

    # print counts
    number_of_virions = sum(counts)

    n = int(number_of_virions * percent)
    np = 0

    for c in counts:
        if n > 0:
            n -= c
            np += 1
    return np


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    fasta = list(SeqIO.parse(args.input, 'fasta'))
    epi = Epistasis()
    epi.build_epistasis_windows_test(fasta)

    print "%s\t%i\t%i\t%i\t%i" % (args.input.name, len(epi.windows),
                                  get_number_of_variants_percentile(fasta),
                                  get_number_of_variants_percentile(fasta, 0.75),
                                  get_number_of_variants_percentile(fasta, 0.9))