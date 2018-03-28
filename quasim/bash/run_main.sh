#!/usr/bin/env bash

export PYTHONPATH='$PYTHONPATH:/home/code/PyBioTools'

. activate smrtenv

OUT=${3}${5}
python /home/code/PyBioTools/quasim/main.py -i $1 -T $2 -L 1000 -o ${OUT} -n ${4}

mkdir ${OUT}/subsamples

for i in ${OUT}/samples/*
do
  subsample_counter.py -r ${i} > ${i%.*}_s.fas -n 50
  mv ${i%.*}_s.fas ${OUT}/subsamples/
done

ARGS=''
for i in ${OUT}/samples/*
do
  ARGS=${ARGS}' '${i}
done

merge_consensuses.py $ARGS > ${OUT}/consensuses.fa

zip -r /var/www/html/quasim/T${2}_n${4}_out${5}.zip ${OUT}/*