#!/usr/bin/env bash

export PYTHONPATH='$PYTHONPATH:/home/code/PyBioTools'

. activate smrtenv

python /home/code/PyBioTools/quasim/intra.py -i $1 -T 500 -L 5000 -o $2


