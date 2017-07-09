#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 1/15/2017
#
######################

import networkx as nx
import argparse
import sys

def __rec_process_root(G, reduced_G, root):
    for s in G.neighbors_iter(root):
        ss = s
        w = 1.0
        while G.out_degree(ss) == 1:
            ss = G.neighbors(ss)[0]
            w += 1
        reduced_G.add_edge(root, ss, weight=w)
        if G.out_degree(ss) == 0:
            continue
        __rec_process_root(G, reduced_G, ss)

def reduce_graph(G):
    root = __find_root(G)

    rg = nx.DiGraph()

    __rec_process_root(G, rg, root)

    return rg

def __find_root(G):
    for n in G:
        if G.in_degree(n) == 0:
            return n
    return None


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-o", dest="output", type=argparse.FileType('w+'), default=sys.stdout)

    args = parser.parse_args()

    G = nx.read_pajek(args.input)

    rg = reduce_graph(G)

    nx.write_pajek(rg, args.output)