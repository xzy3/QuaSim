#!/usr/bin/env bash


PREFIX=$1
OUTDIR=$2

mkdir -p ${OUTDIR}

for i in {1..100}
do
    P1=${PREFIX}/${i}/1/0_t1000.fas
    P2=${PREFIX}/${i}/2/0_t1000.fas
    if [ -f "${P1}" ] && [ -f "${P2}" ]; then
        mkdir -p ${OUTDIR}/${i}
        cp ${P1} ${OUTDIR}/${i}/patient1.fas
        cp ${P2} ${OUTDIR}/${i}/patient2.fas
    fi
done