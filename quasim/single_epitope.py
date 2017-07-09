#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 9/15/2015
#
######################

from Bio import (SeqIO, SeqRecord, Seq)
from Levenshtein import distance

import itertools
import networkx as nx
import argparse
import sys
import os
from quasim import *

# LIST OF DEFAULT PARAMETERS INTRA
T_SIMULATION_TIME=100
D_IMMUNE_SYSTEM_DELAY=5
L_LIVER_SIZE=5000
CDR_INFECTED_CELL_DEATH_RATE=0.1
CIR_PROBABILITY_TO_INFECT_LIVER_CELL=0.1
BCR_B_CELL_RATE=5E-4
DEFAULT_REGION_LENGTH=60

###########
MUTATION_RATE=1E-5
MUTATE_INITIAL_SEQUENCE=0
###########


# LIST OF DEFAULT PARAMETERS INTER
EXP_FOR_SCALE_FREE_GRAPH=2.5

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

def mutate_initial(initial, n_mut=0):
    if n_mut == 0: return initial

    l = len(initial.seq)

    m_pos = random.sample(range(l), n_mut)

    mseq = list(str(initial.seq))

    for m in m_pos:
        mseq[m] = random.choice(Variant.gen_map[mseq[m]])

    assert sum(x != y for (x, y) in zip(mseq, initial.seq)) == n_mut
    initial.seq = Seq.Seq(''.join(mseq), Seq.Alphabet.SingleLetterAlphabet())

def simulate_intrahost(ihost, initial):
    vv = Variant(str(initial))
    ihost.blood.infect(Virion(vv))
    ihost.blood.flow()
    ihost.blood.infect(ihost.liver.produce_virions(mr))

def compute_scale_freeness(G):
    assert isinstance(G, nx.Graph)
    n = G.number_of_nodes()

    sf = sum(G.degree(e[0])*G.degree(e[1]) for e in G.edges_iter())

    return float(sf) / ((n-1)**2)

def filter_dead_ancestries(G, vars):
    assert isinstance(G, nx.Graph)

    if len(G) < 10: return G

    alive_vars = filter(lambda v: v.count() > 0, vars)

    G_nodes = set()

    for v in alive_vars:
        if v in G_nodes: continue
        ancestries = [v]
        while len(ancestries) > 0:
            vv = ancestries[0]

            ancestries.remove(vv)
            if vv in G_nodes: continue
            ancestries.extend(G.predecessors(vv))
            G_nodes.add(vv)

    return G.subgraph(G_nodes)

def get_diversity(variants):
    variants = filter(lambda x: x.count() > 0, variants)
    N = sum(v.count() for v in variants)
    L = max(len(v.seq) for v in variants)
    if len(variants) <= 1:
        return 0.0
    div = sum(s.count() * t.count() * distance(str(s.seq), str(t.seq)) for (s, t) in itertools.combinations(variants, 2))
    div = 2.0 * div / ( N * N * L )
    return div

def get_divergence(variants, initial):
    variants = filter(lambda x: x.count() > 0, variants)
    N = sum(v.count() for v in variants)
    L = max(len(v.seq) for v in variants)
    if len(variants) <= 1:
        return 0.0
    div = sum(s.count() * distance(str(s.seq), str(initial.seq)) for s in variants)
    div = float(div) / (N * L)
    return div

###################
# MAIN METHOD     #
###################
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-T", dest='T', type=int, default=T_SIMULATION_TIME)
    parser.add_argument("-d", dest='delay', type=int, default=D_IMMUNE_SYSTEM_DELAY)
    parser.add_argument("-L", dest='liver_size', type=int, default=L_LIVER_SIZE)
    parser.add_argument("-mr", dest="mutation_rate", type=float, default=MUTATION_RATE)
    parser.add_argument("-cdr", dest='cell_death_rate', type=float, default=CDR_INFECTED_CELL_DEATH_RATE)
    parser.add_argument("-bcr", dest="b_cell_rate", type=float, default=BCR_B_CELL_RATE)
    parser.add_argument("-min_ab", dest="min_ab_eff", type=float, default=None)
    parser.add_argument("-mut", dest="mutate_initial", type=int, default=MUTATE_INITIAL_SEQUENCE)
    parser.add_argument("-o", dest='out_dir', type=str, required=True)
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    cir = CIR_PROBABILITY_TO_INFECT_LIVER_CELL

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

    fasta = list(SeqIO.parse(args.input, 'fasta'))

    initial = random.choice(fasta)

    mutate_initial(initial, args.mutate_initial)
    Virion.epitopes = Epitopes(len(initial.seq), DEFAULT_REGION_LENGTH)

    S = [0]
    hosts = {x: Host(liver_size, initial, args.min_ab_eff) for x in S}
    I = []
    i = random.choice(S)
    ihost = hosts[i]

    vv = Variant(str(initial.seq), ihost)
    ihost.blood.infect(Virion(vv))
    S.remove(i)
    I.append(i)

    print "Initial node ID: %i" % i

    log = open(out_dir_f % "stats.csv", mode="w+")
    log.write("tick\t"
              "diversity\t"
              "divergence\t"
              "# of virions\t"
              "# of variants\n"
              )

    for j in xrange(T):
        tI = []
        for i in I:
            log.write("%i\t" % j)
            ihost = hosts[i]
            #########################
            ihost.tick(j, mr, cdr, bcr, delay, cir)
            v_count = sum(i.count() for i in ihost.blood.variants)
            variants = sum(1 if i.count() > 0 else 0 for i in ihost.blood.variants)

            G = ihost.liver.gen_tree

            print "%i : %i : %i : %.2f%%" % (j, v_count, variants, ihost.liver.infected_portion())
            if v_count == 0:
                print "Population died out"
                sys.exit(-1)
            else:
                log.write("%.3e\t" % get_diversity(ihost.blood.variants))
                log.write("%.3e\t" % get_divergence(ihost.blood.variants, ihost.initial))
                log.write("%i\t" % v_count)
                log.write("%i\t" % variants)
                # log.write("%.3e\t" % compute_scale_freeness(G))
                # log.write("%.3e\t" % compute_scale_freeness(filter_dead_ancestries(G, ihost.blood.variants)))
                log.write('\n')
        # if (j + 1) % dT == 0:
        print_samples(I, hosts, j+1)

        I.extend(tI)
    nx.write_pajek(ihost.liver.gen_tree, out_dir_f % "gen_tree.net")
    nx.write_edgelist(ihost.liver.gen_tree, out_dir_f % "gen_tree.edges")
    nx.write_pajek(filter_dead_ancestries(ihost.liver.gen_tree, ihost.blood.variants), out_dir_f % "gen_tree_alive.net")


