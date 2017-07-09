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
import os
from quasim import *

# LIST OF DEFAULT PARAMETERS INTRA
from quasim.epistasis import Epistasis

T_SIMULATION_TIME=100
D_IMMUNE_SYSTEM_DELAY=5
L_LIVER_SIZE=5000
CDR_INFECTED_CELL_DEATH_RATE=0.1
CIR_PROBABILITY_TO_INFECT_LIVER_CELL=0.1
BCR_B_CELL_RATE=5E-4
DEFAULT_REGION_LENGTH=60

###########
MUTATION_RATE=1E-4
###########


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
    parser.add_argument("-cdr", dest='cell_death_rate', type=float, default=CDR_INFECTED_CELL_DEATH_RATE)
    parser.add_argument("-bcr", dest="b_cell_rate", type=float, default=BCR_B_CELL_RATE)
    parser.add_argument("-min_ab", dest="min_ab_eff", type=float, default=None)
    parser.add_argument("-pm", dest="pm", type=float, default=DEFAULT_ANTIBODY_EFFECTIVENESS)
    parser.add_argument("-pmm", dest="pmm", type=float, default=AB_EFFECTIVENESS_NEIGHBOR)
    parser.add_argument("-o", dest='out_dir', type=str, required=True)
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
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

    # Read input and prepare simulation
    fasta = list(SeqIO.parse(args.input, 'fasta'))

    Variant.epistasis = Epistasis(fasta)

    initial = random.choice(fasta)

    Virion.epitopes = Epitopes(len(initial.seq), DEFAULT_REGION_LENGTH)

    ihost = Host(liver_size, initial, args.min_ab_eff)

    vv = Variant(str(initial.seq), ihost)
    ihost.blood.infect(Virion(vv))


    log = sys.stderr
    log.write("T\t"
              "# of virions\t"
              "# of variants\t")

    for j in xrange(T):
        tI = []

        #########################
        ihost.tick(j, mr, cdr, bcr, delay, cir, pm=pm, pmm=pmm)
        v_count = sum(i.count() for i in ihost.blood.variants)
        variants = sum(1 if i.count() > 0 else 0 for i in ihost.blood.variants)

        G = ihost.liver.gen_tree

        if v_count == 0:
            print "Population died out"
            sys.exit(-1)
        else:
            log.write("%i : %i : %i : %.2f%%" % (j, v_count, variants, ihost.liver.infected_portion()))
            log.write('\n')
        print_sample(ihost, j+1)




