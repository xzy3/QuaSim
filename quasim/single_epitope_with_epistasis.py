#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 3/17/2016
#
######################

from Bio import (SeqIO, SeqRecord, Seq)

import argparse
import sys
import networkx as nx
import os
from quasim import *

# LIST OF DEFAULT PARAMETERS INTRA
from quasim.epistasis import Epistasis
from quasim.single_epitope_with_profile import (get_diversity, get_divergence)


def print_sample(ihost, j, count_threshold=0):
    i = 0
    variants = filter(lambda v: v.count() > count_threshold, ihost.blood.variants)
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


###################
# MAIN METHOD     #
###################
if __name__=='__main__':

    # Setup arguments and parse them

    parser = argparse.ArgumentParser()
    parser.add_argument("-T", dest='T', type=int, default=T_SIMULATION_TIME)
    parser.add_argument("-d", dest='delay', type=int, default=D_IMMUNE_SYSTEM_DELAY)
    parser.add_argument("-L", dest='liver_size', type=int, default=L_LIVER_SIZE)
    parser.add_argument("-mr", dest="mutation_rate", type=float, default=MUTATION_RATE)
    parser.add_argument("-S", dest='start', type=int, default=0)
    parser.add_argument("-cdr", dest='cell_death_rate', type=float, default=CDR_INFECTED_CELL_DEATH_RATE)
    parser.add_argument("-bcr", dest="b_cell_rate", type=float, default=BCR_B_CELL_RATE)
    parser.add_argument("-min_ab", dest="min_ab_eff", type=float, default=None)
    parser.add_argument("-pm", dest="pm", type=float, default=DEFAULT_ANTIBODY_EFFECTIVENESS)
    parser.add_argument("-pmm", dest="pmm", type=float, default=AB_EFFECTIVENESS_NEIGHBOR)
    parser.add_argument("-o", dest='out_dir', type=str, required=True)
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument("--input-seqs", dest='input_seqs', type=argparse.FileType('r'), default=None)
    args = parser.parse_args()

    cir = CIR_PROBABILITY_TO_INFECT_LIVER_CELL

    # Prepare output directory

    if os.path.isdir(args.out_dir):
        os.system('rm -r %s' % args.out_dir)
    if os.path.exists(args.out_dir):
        os.remove(args.out_dir)
    os.makedirs(args.out_dir)
    out_dir_f = args.out_dir + '/%s'

    seed = None

    # Transfer input parameters
    num = 1
    liver_size = args.liver_size
    delay = args.delay

    mr = args.mutation_rate

    # cp = 0.02
    cdr = args.cell_death_rate
    bcr = args.b_cell_rate
    Variant.creact = True
    T = args.T

    pm = args.pm
    pmm = args.pmm

    start = args.start

    # Read input and prepare simulation
    fasta = list(SeqIO.parse(args.input, FASTA))
    input_seqs = fasta if args.input_seqs is None else list(SeqIO.parse(args.input_seqs, FASTA))

    Variant.epistasis = Epistasis(fasta + input_seqs)

    initial = random.choice(input_seqs)

    Virion.epitopes = Epitopes(len(initial.seq), DEFAULT_REGION_LENGTH)

    ihost = Host(liver_size, initial, args.min_ab_eff)

    vv = Variant(str(initial.seq), ihost)
    ihost.blood.infect(Virion(vv))

    log = open(out_dir_f % "stats.csv", mode="w+")
    log.write("tick\t"
              "diversity\t"
              "divergence\t"
              "# of virions\t"
              "# of variants\n"
              )

    for j in xrange(start, T):
        tI = []

        #########################
        ihost.tick(j, mr, cdr, bcr, delay, cir, pm=pm, pmm=pmm)
        v_count = sum(i.count() for i in ihost.blood.variants)
        variants = sum(1 if i.count() > 0 else 0 for i in ihost.blood.variants)

        G = ihost.liver.gen_tree

        print "%i : %i : %i : %.2f%%" % (j, v_count, variants, ihost.liver.infected_portion())
        if v_count == 0:
            sys.stderr.write("Population died out\n")
            sys.exit(-1)
        else:
            log.write("%.3e\t" % get_diversity(ihost.blood.variants))
            log.write("%.3e\t" % get_divergence(ihost.blood.variants, ihost.initial))
            log.write("%i\t" % v_count)
            log.write("%i\t" % variants)
            # log.write("%.3e\t" % compute_scale_freeness(G))
            # log.write("%.3e\t" % compute_scale_freeness(filter_dead_ancestries(G, ihost.blood.variants)))
            log.write('\n')
        print_sample(ihost, j+1)




