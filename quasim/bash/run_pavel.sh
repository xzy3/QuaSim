#!/usr/bin/env bash

for i in {1..10}
do
   /home/aartyomenko/sasha/bin/python /home/aartyomenko/BioTools/quasim/main.py -o /data1/Sasha/Pavel/${i} -i ~/BioTools/quasim/initial.fas -n 20 -T 300
done

for i in {11..20}
do
   /home/aartyomenko/sasha/bin/python /home/aartyomenko/BioTools/quasim/main.py -o /data1/Sasha/Pavel/${i} -i ~/BioTools/quasim/initial.fas -n 45 -T 300
done

for i in {21..30}
do
   /home/aartyomenko/sasha/bin/python /home/aartyomenko/BioTools/quasim/main.py -o /data1/Sasha/Pavel/${i} -i ~/BioTools/quasim/initial.fas -n 70 -T 400
done

for i in {31..40}
do
   /home/aartyomenko/sasha/bin/python /home/aartyomenko/BioTools/quasim/main.py -o /data1/Sasha/Pavel/${i} -i ~/BioTools/quasim/initial.fas -n 95 -T 400
done

for i in {41..50}
do
   /home/aartyomenko/sasha/bin/python /home/aartyomenko/BioTools/quasim/main.py -o /data1/Sasha/Pavel/${i} -i ~/BioTools/quasim/initial.fas -n 120 -T 500
done
