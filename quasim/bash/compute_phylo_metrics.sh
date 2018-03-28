#!/usr/bin/env bash

PREFIX="/research_data/sasha/CDC/454_filtered/Acutes_NGS"
OUTLOG="/research_data/sasha/CDC/454_filtered/NGS_acute.log"


#PREFIX="/research_data/sasha/CDC/454_filtered/Chronics_NGS"
#OUTLOG="/research_data/sasha/CDC/454_filtered/NGS_chronic.log"

export PYTHONPATH=$PYTHONPATH:/home/code/PyBioTools

. activate smrtenv

printf "" > ${OUTLOG}

for i in ${PREFIX}/*.fas
    do
    if [ -f "${i%.*}.nwk" ]; then
        /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/phylo.py -i ${i%.*}.nwk >> ${OUTLOG}
    fi
    if [ ! -f "${i%.*}.nwk" ]; then
        /home/code/anaconda/envs/smrtenv/bin/python -u /home/code/PyBioTools/quasim/phylo.py -i ${i} -t ${i%.*}.nwk >> ${OUTLOG}
    fi

    echo ${i}
    done

echo "done!"