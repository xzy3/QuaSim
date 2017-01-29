#! /usr/bin/env python

from Bio import (SeqIO, SeqRecord, Seq)
from Levenshtein import distance

import itertools
import argparse
import sys
import os
from quasim import *


def print_samples(I, hosts, j):
    for i in I:
        ihost = hosts[i]
        variants = filter(lambda v: v.count() > 0, ihost.blood.variants)
        seqs = map(lambda (c, v):
                   SeqRecord.SeqRecord(Seq.Seq(v.seq, Seq.Alphabet.SingleLetterAlphabet()),
                                       id="%i_%i_%i" % (j, c, v.count()), description=''),
                   enumerate(sorted(variants, key=lambda x: -x.count())))
        if len(seqs) > 0:
            f = open(out_dir_f % "%i_t%i.fas" % (i, j), 'w+')
            SeqIO.write(seqs, f, 'fasta')

def simulate_intrahost(ihost, initial):
    vv = Variant(str(initial))
    ihost.blood.infect(Virion(vv))
    ihost.blood.flow()
    ihost.blood.infect(ihost.liver.produce_virions(mr))

def get_diversity(variants):
    div = 0.0
    variants = filter(lambda x: x.count() > 0, variants)
    N = sum(v.count() for v in variants)
    L = max(len(v.seq) for v in variants)
    if len(variants) <= 1:
        return 0.0
    for (s, t) in itertools.combinations(variants, 2):
        div += s.count() * t.count() * distance(str(s.seq), str(t.seq))
    div = 2 * div / ( N * N * L )
    return div

def get_divergence(variants, initial):
    div = 0.0
    variants = filter(lambda x: x.count() > 0, variants)
    N = sum(v.count() for v in variants)
    L = max(len(v.seq) for v in variants)
    if len(variants) <= 1:
        return 0.0
    for s in variants:
        div += s.count() * distance(str(s.seq), str(initial.seq))
    div = div / (N * L)
    return div

###################
# MAIN METHOD     #
###################
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-T", dest='T', type=int, default=20)
    parser.add_argument("-d", dest='delay', type=int, default=5)
    parser.add_argument("-L", dest='liver_size', type=int, default=1000)
    parser.add_argument("-mr", dest="mutation_rate", type=float, default=7.5E-6)
    parser.add_argument("-cdr", dest='cell_death_rate', type=float, default=0.01)
    parser.add_argument("-bcr", dest="b_cell_rate", type=float, default=2E-4)
    parser.add_argument("-o", dest='out_dir', type=str, required=True)
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    cir = 0.1

    if os.path.isdir(args.out_dir):
        os.system('rm -r %s' % args.out_dir)
    if os.path.exists(args.out_dir):
        os.remove(args.out_dir)
    os.makedirs(args.out_dir)
    out_dir_f = args.out_dir + '/%s'


    seed = None

    # Input parameters
    num = 1
    liver_size = args.liver_size   #5000
    delay = args.delay #10
    exp = 2.5

    mr = args.mutation_rate #0.005

    # cp = 0.02
    cdr = args.cell_death_rate    #0.25
    bcr = args.b_cell_rate    #0.05
    Variant.creact = True
    T = args.T
    dT = 20 # step in cycles before printing output

    fasta = SeqIO.parse(args.input, 'fasta')

    initial = fasta.next()

    print initial
    Virion.epitopes = Epitopes(len(initial.seq), 60)

    S = [0]
    hosts = { x: Host(liver_size, initial) for x in S}
    I = []
    i = random.choice(S)
    ihost = hosts[i]

    vv = Variant(str(initial.seq), ihost)
    ihost.blood.infect(Virion(vv))
    S.remove(i)
    I.append(i)

    print "Initial node ID: %i" % i

    log = open(out_dir_f % "stats.csv", mode="w+")
    log.write("tick\tdiversity\tdivergence\t# of virions\t# of variants\n")

    for j in xrange(T):
        tI = []
        for i in I:
            log.write("%i\t" % j)
            ihost = hosts[i]
            #########################
            ihost.tick(mr, cdr, bcr, delay, cir)
            v_count = sum(i.count() for i in ihost.blood.variants)
            variants = sum(1 if i.count() > 0 else 0 for i in ihost.blood.variants)

            print "%i : %i" % (j, v_count)
            if v_count == 0:
                print "Population died out"
                sys.exit(-1)
            log.write("%.3e\t" % get_diversity(ihost.blood.variants))
            log.write("%.3e\t" % get_divergence(ihost.blood.variants, ihost.initial))
            log.write("%i\t" % v_count)
            log.write("%i\n" % variants)
        # if (j + 1) % dT == 0:
        print_samples(I, hosts, j+1)

        I.extend(tI)
