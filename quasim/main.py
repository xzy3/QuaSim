#! /usr/bin/env python
import networkx as nx
import networkx.drawing.nx_pydot as nd
from networkx.utils import (powerlaw_sequence, create_degree_sequence)
from Bio import (SeqIO, SeqRecord, Seq)

import argparse
import sys
import os
from quasim import *


def print_samples(I, hosts):
    for i in I:
        ihost = hosts[i]

        seqs = map(lambda (c, v):
                   SeqRecord.SeqRecord(Seq.Seq(v.seq, Seq.Alphabet.SingleLetterAlphabet()),
                                       id="%i_%i" % (i, c), description="%i" % v.count()),
                   enumerate(sorted(ihost.blood.variants, key=lambda x: -x.count())))
        if len(seqs) > 0:
            f = open(out_dir_f % "samples/%i.fas" % i, 'w+')
            SeqIO.write(seqs, f, 'fasta')

###################
# MAIN METHOD     #
###################
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", dest='num', type=int, default=100)
    parser.add_argument("-T", dest='T', type=int, default=10)
    parser.add_argument("-d", dest='delay', type=int, default=2)
    parser.add_argument("-L", dest='liver_size', type=int, default=1000)
    parser.add_argument("-mr", dest="mutation_rate", type=float, default=1.5E-3)
    parser.add_argument("-cdr", dest='cell_death_rate', type=float, default=0.025)
    parser.add_argument("-bcr", dest="b_cell_rate", type=float, default=1E-4)
    parser.add_argument("-o", dest='out_dir', type=str, required=True)
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    if os.path.isdir(args.out_dir):
        os.system('rm -r %s' % args.out_dir)
    if os.path.exists(args.out_dir):
        os.remove(args.out_dir)
    os.makedirs(args.out_dir)
    out_dir_f = args.out_dir + '/%s'


    seed = None

    # Input parameters
    num = args.num
    liver_size = args.liver_size   #5000
    delay = args.delay #10
    exp = 2.5

    mr = args.mutation_rate #0.005

    cp = 0.005
    cdr = args.cell_death_rate    #0.25
    bcr = args.b_cell_rate    #0.05
    Variant.creact = True
    T = args.T

    fasta = SeqIO.parse(args.input, 'fasta')

    initial = fasta.next()

    print initial

    # Generate network
    # noinspection PyDeprecation
    sequence = create_degree_sequence(num, powerlaw_sequence, exponent=exp)
    graph = nx.configuration_model(sequence, seed=seed)
    loops = graph.selfloop_edges()
    # remove parallel edges and self-loops
    graph = nx.Graph(graph)
    graph.remove_edges_from(loops)
    # get largest connected component
    # unfortunately, the iterator over the components is not guaranteed to be sorted by size
    component = max(nx.connected_components(graph), key=len)
    lcc = graph.subgraph(component)

    nd.write_dot(lcc, out_dir_f % "scale_free.dot")

    print "Number of nodes in generated scale free graph: %i" % len(lcc)

    transmission = nx.DiGraph()
    # transmission.add_nodes_from(lcc.nodes())
    Virion.epitopes = Epitopes(len(initial.seq), k=5)
    S = lcc.nodes()
    hosts = { x: Host(liver_size) for x in S}
    I = []
    i = random.choice(S)
    ihost = hosts[i]

    vv = Variant(str(initial.seq), ihost)
    ihost.blood.infect(Virion(vv))
    S.remove(i)
    I.append(i)

    print "Initial node ID: %i" % i

    log = open(out_dir_f % "time.log", mode="w+")
    log_format = "%i - %i\n"
    log.write(log_format % (i, 0))

    for j in xrange(T):
        tI = []
        if T-j < 10: cp = 0.0
        for i in I:
            ihost = hosts[i]
            #########################
            ihost.tick(mr, cdr, bcr, delay)
            v_count = sum(i.count() for i in ihost.blood.variants)

            print "%i : %i" % (j, v_count)
            if v_count == 0:
                print "Population died out"
                sys.exit(-1)
            for ti in nx.all_neighbors(lcc, i):
                if ihost.rnd.random() > cp: continue
                if ti not in S: continue
                v = ihost.rnd.choice(tuple(ihost.blood.variants))
                S.remove(ti)
                target = hosts[ti]

                var = Variant(v.seq, target)
                virion = Virion(var)

                target.blood.infect(virion)
                log.write(log_format % (ti, j))
                tI.append(ti)
                transmission.add_edge(i, ti)
                break

        I.extend(tI)

    nd.write_dot(transmission, out_dir_f % "transmission.dot")
    os.makedirs(out_dir_f % 'samples')

    print_samples(I, hosts)

