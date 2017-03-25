#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 3/15/2017
#
######################

import argparse
import sys
import os
import re
import random

from Bio import (AlignIO, SeqRecord, Seq)
from Bio.Align import MultipleSeqAlignment

from quasim import phylo
from collections import Counter
########################################################

pattern="[0-9]*$"
BOOTSTRAP_DEFAULT_SIZE=40
TREE_SIZE_DEFAULT_THRESHOLD=20

def compute_p_value(seqs, N, n, function, out=sys.stderr):
    """
    Compute empirical p_value using n bootstraps
    N is a sample size
    function returns computed property e. g. Sackin index
    """
    return compute_p_values(seqs, N, n, [function], out)[0]

def compute_p_values(seqs, N, n, functions, out=sys.stderr):
    """
    Compute empirical p_values using n bootstraps
    N is a sample size
    functions (list) each returns computed property e. g. Sackin index
    """

    seqs_with_counts = Counter({x: get_count(seqs[x]) for x in seqs})
    seqs_array = list(seqs_with_counts.elements())

    real_values = [function(real_tree) for function in functions]

    n_less = [0.0 for _ in functions]
    n_greater = [0.0 for _ in functions]

    def generate_random_seq(seqs, seqs_array, l, id):
        str_seq = "".join(seqs[random.choice(seqs_array)].seq[i] for i in range(l))
        return SeqRecord.SeqRecord(Seq.Seq(str_seq, alphabet=Seq.Alphabet.SingleLetterAlphabet()), id=str(id))

    l = len(seqs[seqs_array[0]])

    for i in range(n):
        selected_sim_seqs = set()
        id = 0
        while len(selected_sim_seqs) < N:
            selected_sim_seqs.add(generate_random_seq(seqs, seqs_array, l, id))
            id += 1

        msa_sim_seqs = MultipleSeqAlignment([x for x in selected_sim_seqs])
        sim_tree = phylo.build_phylogenetic_tree(msa_sim_seqs)
        try:
            sim_tree.root_at_midpoint()
        except Exception:
            continue
        boot_values = [function(sim_tree) for function in functions]

        boot_vs_real_values = zip(boot_values, real_values)

        boots_less = (boot_value <= real_value for boot_value, real_value in boot_vs_real_values)
        n_less = [sum(x) for x in zip(n_less, boots_less)]

        boots_greater = (boot_value >= real_value for boot_value, real_value in boot_vs_real_values)
        n_greater = [sum(x) for x in zip(n_greater, boots_greater)]

        out.write("%i\t" % (i+1))
        out.write("\t".join("%i %i" % (real_value, boot_value) for boot_value, real_value in boot_vs_real_values))
        out.write("\n")

    return [min(2 * min(n_l, n_g) / n, 1.0) for n_l, n_g in zip(n_less, n_greater)]

def get_count(fas):
    m = re.search(pattern, fas.description)
    if m is not None:
        return int(m.group(0))
    return 1

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input_real', type=argparse.FileType('r'), required=True)
    parser.add_argument("-s", dest='input_sim', type=argparse.FileType('r'), required=True)
    parser.add_argument("-n", dest='boot', type=int, default=BOOTSTRAP_DEFAULT_SIZE)
    parser.add_argument("-t", dest='min_tree_size', type=int, default=TREE_SIZE_DEFAULT_THRESHOLD)
    parser.add_argument("-p", dest='pattern', type=str, default=pattern)
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    in_sequences = {x.id: x for x in AlignIO.read(args.input_real, phylo.FASTA)}
    sim_in_sequences = {x.id: x for x in AlignIO.read(args.input_sim, phylo.FASTA)}
    real_tree = phylo.build_phylogenetic_tree(MultipleSeqAlignment(in_sequences.values()))

    N = len(in_sequences)

    if N < args.min_tree_size:
        sys.stderr.write("Tree size is too small!\n")
        sys.exit(-1)

    try:
        real_tree.root_at_midpoint()
    except:
        sys.stderr.write("Error processing tree!\n")
        sys.exit(-1)

    p_functions = [phylo.get_sackin_index, phylo.get_colless_index, phylo.count_cherries, phylo.get_shape_hop_length]
    args.output.write("%s\t" % os.path.basename(args.input_real.name))
    args.output.write("%s\t" % "\t".join("%.3f" % p for p in compute_p_values(sim_in_sequences, N, args.boot, p_functions)))
    args.output.write("%i\n" % N)