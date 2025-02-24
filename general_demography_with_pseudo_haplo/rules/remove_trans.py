#!/bin/python3

import sys
import csv

IN_FREQ_FILE=snakemake.input[0]
OUT_FILE=snakemake.output[0]
out_no_trans=csv.writer(open(OUT_FILE,'w'), delimiter='\t')
transversion=[{'A','C'},{'A','T'},{'G','C'},{'G','T'}]

with open(IN_FREQ_FILE) as tsv:
	for line in csv.reader(tsv, dialect="excel", delimiter=' ', skipinitialspace=True):
#CHANGED THIS LINE CAUSE WAS GIVING ERROR. DONT KNOW WHY BUT MY FREQ FILE IS SEPARATED BY MULTIPLE SPACES
		alleles=set(line[2]+line[3])
		if alleles in transversion:
			out_no_trans.writerow(line)
