#!/usr/bin/env bash

export PYTHONPATH=$PYTHONPATH:/home/code/PyBioTools

OUTNAME=$2
PREFIXNAME=$3
PREFIX="/research_data/sasha/CDC/EPLD/"${PREFIXNAME}
OUTLOG="/research_data/sasha/CDC/EPLD/"${OUTNAME}
SIM_PATH=$1

if [ ! -f ${SIM_PATH} ]; then
    echo "Path: "${SIM_PATH}" is incorrect!"
    exit 1
fi

. activate smrtenv

printf "" > ${OUTLOG}

for i in ${PREFIX}/*.fas
do
#    printf $(basename ${i})"\t" >> ${OUTLOG}
    if [ -f "${i%.*}.nwk" ]; then
        /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/significance_level.py -i ${i} -s ${SIM_PATH} >> ${OUTLOG}
    fi
    if [ ! -f "${i%.*}.nwk" ]; then
        /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/significance_level.py -i ${i%.*}.fas -s ${SIM_PATH} >> ${OUTLOG}
    fi

    echo ${i}
done