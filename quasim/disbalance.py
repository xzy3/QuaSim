#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 1/11/2017
#
######################

import argparse
import re

from Bio import (SeqIO)

from quasim.entropy import Profile


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

    fasta = list(SeqIO.parse(args.input, 'fasta'))

    k = args.k
    patt = args.pattern

    profile = Profile()
    profile.load_from_fasta(fasta, k, patt)

    m, std = profile.std()
    print '%s\t%f\t%f\t%i\t%i' % (args.input.name, m, std, profile.size, len(profile.seqs))
    # print '\n'.join(map(str, profile.get_positional_entropy()))