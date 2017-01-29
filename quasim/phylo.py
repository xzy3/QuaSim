#! /usr/bin/env python

from Bio.Phylo.TreeConstruction import (DistanceTreeConstructor, DistanceCalculator)
from Bio import (AlignIO, Phylo)
import networkx.drawing.nx_pydot as nd
import pylab

if __name__=='__main__':
    aln = AlignIO.read('/home/code/PyBioTools/quasim/data/tsts2/subsamples/9_s.fas', 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    print dm
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    net = Phylo.to_networkx(tree)
    nd.write_dot(net, 'test.dot')

