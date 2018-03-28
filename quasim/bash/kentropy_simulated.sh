#!/usr/bin/env bash

PREF=$1

. activate smrtenv

export PYTHONPATH=$PYTHONPATH:/home/code/PyBioTools

for i in {50..1000..1}
do
    INPUT=${PREF}/0_t${i}.fas
    if [ -f ${INPUT} ]; then
        printf ${i}
        /home/code/anaconda2/envs/smrtenv/bin/python2.7 -u /home/code/PyBioTools/quasim/kentropy_table.py -i ${INPUT}
    fi
done > ${PREF}/log.log
