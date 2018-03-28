#!/usr/bin/env bash

I=$1
ST=500
RT=200
TRT=300

export PYTHONPATH=$PYTHONPATH:/home/code/PyBioTools

/home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/single_epitope.py -i /home/code/PyBioTools/quasim/data/initial.fas -o /home/code/PyBioTools/quasim/data/sergey/${I} -T ${ST} -L 1000

/home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/single_epitope.py -i /home/code/PyBioTools/quasim/data/sergey/${I}/0_t${TRT}.fas -o /home/code/PyBioTools/quasim/data/sergey/${I}/sub -T ${RT} -L 1000

if [ -f "/home/code/PyBioTools/quasim/data/sergey/${I}/0_t${ST}.fas" ] && [ -f "/home/code/PyBioTools/quasim/data/sergey/${I}/sub/0_t${RT}.fas" ]
then
    mkdir -p /home/code/PyBioTools/quasim/data/sergey/res/${I}
    cp /home/code/PyBioTools/quasim/data/sergey/${I}/0_t${ST}.fas /home/code/PyBioTools/quasim/data/sergey/res/${I}/source.fas
    cp /home/code/PyBioTools/quasim/data/sergey/${I}/sub/0_t${RT}.fas /home/code/PyBioTools/quasim/data/sergey/res/${I}/recipient.fas
fi
