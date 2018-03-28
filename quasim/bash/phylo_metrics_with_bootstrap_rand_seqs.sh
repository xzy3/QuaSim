#!/usr/bin/env bash

export PYTHONPATH=$PYTHONPATH:/home/code/PyBioTools

OUTNAME=$1
PREFIXNAME=$2
SIM_PATH=$3
PREFIX="/research_data/sasha/CDC/EPLD/"${PREFIXNAME}
OUTLOG="/research_data/sasha/CDC/EPLD/"${OUTNAME}


. activate smrtenv

printf "" > ${OUTLOG}

for i in ${PREFIX}/*.fas
do
#    printf $(basename ${i})"\t" >> ${OUTLOG}
    if [ -f "${i%.*}.nwk" ]; then
        /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/random_sequences.py -i ${i} -s ${SIM_PATH} >> ${OUTLOG}
    fi
    if [ ! -f "${i%.*}.nwk" ]; then
        /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/random_sequences.py -i ${i%.*}.fas -s ${SIM_PATH} >> ${OUTLOG}
    fi

    echo ${i}
done