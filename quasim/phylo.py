#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 1/15/2017
#
######################

#! /usr/bin/env python

from Bio.Phylo.TreeConstruction import (DistanceTreeConstructor, DistanceCalculator)
from Bio import (AlignIO, Phylo)
from Bio.Phylo import NewickIO
from quasim.disbalance import get_count
import argparse
import sys

DISTANCE_TYPE='identity'
TREE_CONSTRUCTION_ALGORITHM='nj'
FASTA='fasta'
FASTA_EXTENSIONS=['.fasta', '.fas', '.fa']
NEWICK_EXTENSIONS=['.nwk', '.newick']
MIN_COUNT=1
NUM_OF_VIRIONS=0

def build_phylogenetic_tree(seqs):
    calculator = DistanceCalculator(DISTANCE_TYPE)
    # Print distance matrix for testing
    # distance_matrix = calculator.get_distance(seqs)

    constructor = DistanceTreeConstructor(calculator, TREE_CONSTRUCTION_ALGORITHM)

    tree = constructor.build_tree(seqs)

    return tree

def read_fasta_or_newick_and_return_tree(path, nwk_path = None, patt = None):
    global NUM_OF_VIRIONS
    if any(path.name.endswith(x) for x in FASTA_EXTENSIONS):
        seqs = AlignIO.read(path, FASTA)
        seqs._records = [x for x in seqs if get_count(x, patt) > MIN_COUNT]
        NUM_OF_VIRIONS = int(sum(get_count(x, patt) for x in seqs))

        if len(seqs) <= 2: return None
        tree = build_phylogenetic_tree(seqs)
        if nwk_path is not None and tree is not None:
            NewickIO.write([tree], nwk_path)
    elif any(path.name.endswith(x) for x in NEWICK_EXTENSIONS):
        tree = NewickIO.parse(path).next()

    # Root the tree if necessary
    if not tree.rooted:
        tree.root_at_midpoint()

    return tree


def get_shape_length(tree):
    clades = tree.get_terminals()
    return sum(tree.distance(clade) for clade in clades)

def get_shape_hop_length(tree):
    clades = tree.find_clades() # all clades

    return sum(len(tree.get_path(clade)) for clade in clades)

def get_sackin_index(tree):
    clades = tree.get_terminals()  # only terminals
    isn = sum(len(tree.get_path(clade)) for clade in clades)
    # n = tree.count_terminals()
    # exp_isn = 2 * n * sum(x ** -1 for x in range(2, n+1)) # E(isn)=2n sum_{j=2..n} {1/j}
    return isn

def count_cherries(tree):
    clades = tree.find_clades()

    return sum(1 if sum(y.is_terminal() for y in x) >= 2 else 0 for x in clades)

def get_colless_index(tree):
    clades = tree.find_clades()

    def size(clade):
        if clade.is_terminal():
            return 1
        return sum(size(c) for c in clade)

    return sum(abs(size(clade[0]) - size(clade[1])) if len(clade) == 2 else 0 for clade in clades)

################################################################################################

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-t", dest='out_tree', type=argparse.FileType('w+'), default=None)
    parser.add_argument("-p", dest='patt', type=str, default="[0-9]*$")
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    patt = args.patt
    tree = read_fasta_or_newick_and_return_tree(args.input, args.out_tree, patt)

    if tree is None:
        # args.output.write("%.4e\t%.4e\t%i\n" % (1, 1, 2))
        sys.exit(0)

    shape_len = get_shape_length(tree)
    total_len = tree.total_branch_length()

    hop_len = get_shape_hop_length(tree)
    sackin = get_sackin_index(tree)
    n = tree.count_terminals()

    hop_r = hop_len * 1.0 / (n * (n-1)) # all clades
    cherries =  count_cherries(tree)
    colless = get_colless_index(tree)


    args.output.write("%.4e\t%.4e\t%.4e\t%d\t%d\t%i\t%i\n" %
                      (shape_len/total_len, hop_r, sackin, cherries, colless, tree.count_terminals(), NUM_OF_VIRIONS))
    # print "Shape len: %.4e and total len: %.4e ratio: %.4e" % (shape_len, total_len, shape_len/total_len)
    # print "Hop length: %.4e ratio: %.4e" % (hop_len, hop_len * 1.0 / (n * (n-1)))