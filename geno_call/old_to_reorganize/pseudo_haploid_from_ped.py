#!/bin/python3

# Importing Modules
import sys
import csv
import random
import itertools as it

# Defining input
in_file=snakemake.input[0]
out_file=snakemake.output[0]

with open(out_file, 'w') as o_f:
    with open(in_file,'r') as in_f:
        for line in csv.reader(in_f, dialect="excel", delimiter=' '):
            out_head=line[0:6]
            in_seq=line[6:]
            in_geno=list(zip(*[iter(in_seq)]*2))
            out_hap=[random.choice(el) for el in in_geno]
            out_geno=list(zip(out_hap,out_hap))
            out_seq=list(it.chain.from_iterable(out_geno))
            out_line=" ".join(out_head+out_seq)
            o_f.write(out_line+"\n")
