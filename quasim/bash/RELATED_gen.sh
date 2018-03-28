#!/usr/bin/env bash

export PYTHONPATH='$PYTHONPATH:/home/code/PyBioTools'

. activate smrtenv

I=/home/code/PyBioTools/quasim/data/initial.fas

mkdir -p /home/code/PyBioTools/quasim/data/RELATED

I1=/home/code/PyBioTools/quasim/data/RELATED/${1}/1
I2=/home/code/PyBioTools/quasim/data/RELATED/${1}/2

T=1000
MIN_AB=0.75
L=1000
OUT=/home/code/PyBioTools/quasim/data/RELATED/log.log

mkdir -p /home/code/PyBioTools/quasim/data/RELATED/${1}

if [ ! -f "${I1}/0_t${L}.fas" ]; then
    /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/single_epitope.py -i ${I} -o ${I1} -T ${T} -min_ab ${MIN_AB} -L ${L}
fi
if [ ! -f "${I2}/0_t${L}.fas" ]; then
    /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/single_epitope.py -i ${I} -o ${I2} -T ${T} -min_ab ${MIN_AB} -L ${L}
fi

/home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/meandist.py -i1 ${I1}/0_t${L}.fas -i2 ${I2}/0_t${L}.fas >> ${OUT}