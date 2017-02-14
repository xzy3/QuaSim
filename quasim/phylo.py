#! /usr/bin/env python

from Bio.Phylo.TreeConstruction import (DistanceTreeConstructor, DistanceCalculator)
from Bio import (AlignIO, Phylo)
from Bio.Phylo import NewickIO
import argparse
import sys

DISTANCE_TYPE='identity'
TREE_CONSTRUCTION_ALGORITHM='nj'
FASTA='fasta'

def build_phylogenetic_tree(seqs):
    calculator = DistanceCalculator(DISTANCE_TYPE)
    # Print distance matrix for testing
    # distance_matrix = calculator.get_distance(seqs)

    constructor = DistanceTreeConstructor(calculator, TREE_CONSTRUCTION_ALGORITHM)

    tree = constructor.build_tree(seqs)

    return tree

def read_fasta_and_return_tree(path):
    seqs = AlignIO.read(path, FASTA)
    if len(seqs) <= 2: return None
    tree = build_phylogenetic_tree(seqs)
    return tree


def get_shape_length(tree):
    clades = tree.get_terminals()
    return sum(tree.distance(clade) for clade in clades)

def get_shape_hop_length(tree):
    clades = tree.find_clades() # all clades
    # clades = tree.get_terminals() # only terminals

    return sum(len(tree.get_path(clade)) for clade in clades)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', type=argparse.FileType('r'), required=True)
    parser.add_argument("-t", dest='out_tree', type=argparse.FileType('w+'), default=None)
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    tree = read_fasta_and_return_tree(args.input)
    if args.out_tree is not None and tree is not None:
        NewickIO.write([tree], args.out_tree)

    if tree is None:
        # args.output.write("%.4e\t%.4e\t%i\n" % (1, 1, 2))
        sys.exit(0)

    shape_len = get_shape_length(tree)
    total_len = tree.total_branch_length()

    hop_len = get_shape_hop_length(tree)
    n = tree.count_terminals()

    hop_r = hop_len * 1.0 / (n * (n-1)) # all clades
    # hop_r = hop_len * 2.0 / (n * (n-1)) # only terminals


    args.output.write("%.4e\t%.4e\t%i\n" % (shape_len/total_len, hop_r, tree.count_terminals()))
    # print "Shape len: %.4e and total len: %.4e ratio: %.4e" % (shape_len, total_len, shape_len/total_len)
    # print "Hop length: %.4e ratio: %.4e" % (hop_len, hop_len * 1.0 / (n * (n-1)))