#! /usr/bin/env python

######################
#
#  Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
#  Created: 1/18/2017
#
######################
import argparse
import re
from Bio import SeqIO
import glob
import sys
from collections import defaultdict

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input_dir', type=str, required=True)
    parser.add_argument("-o", dest='output', type=argparse.FileType('w+'), default=sys.stdout)
    args = parser.parse_args()

    time_series = {}
    output = args.output
    gl_set = set()


    for f in glob.glob("%s/*.fas" % args.input_dir):
        fasta = list(SeqIO.parse(f, 'fasta'))
        m = re.search('t[0-9]*', f)
        index = int(m.group(0)[1:])
        dd = defaultdict(lambda: .0)
        for seq in fasta:
            m = re.search("[0-9]*$", seq.description)
            if m is not None:
                dd[str(seq.seq)] = float(m.group(0))
                gl_set.add(str(seq.seq))
        time_series[index] = dd

    for i in time_series:
        dd = time_series[i]
        n = sum(dd.values())
        for k in dd:
            dd[k] /= n

    output.write("time\t%s\n" % "\t".join(map(str, time_series)))

    for seq in gl_set:
        output.write("%s\t%s\n" % (seq, '\t'.join(map(lambda i: str(time_series[i+1][seq]), xrange(len(time_series))))))

