#!/usr/bin/env bash

I=$1

export PYTHONPATH=$PYTHONPATH:/home/code/PyBioTools

A09=0.95
A099=0.99
OUT_LOG='/home/code/PyBioTools/quasim/data/stats2.log'
/home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/single_epitope.py -i /home/code/PyBioTools/quasim/data/initial.fas -o /home/code/PyBioTools/quasim/data/2k_liver/cir${I} -T 1000 -L 2000 -min_ab ${A09}
/home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/single_epitope.py -i /home/code/PyBioTools/quasim/data/initial.fas -o /home/code/PyBioTools/quasim/data/2k_liver/nocir${I} -T 1000 -L 2000 -min_ab ${A099}

#/home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/reduce_graph.py -i data/09model${I}/gen_tree_alive.net -o data/09model${I}/gen_tree_alive_r.net
#/home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/reduce_graph.py -i data/099model${I}/gen_tree_alive.net -o data/099model${I}/gen_tree_alive_r.net
