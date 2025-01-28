#!/usr/bin/env python3

import argparse
import sys
import csv

parser = argparse.ArgumentParser(prog="count_pairwise_vcftab.py", description="A program to count pairwise variant differences between strains")
parser.add_argument("--input",help="The input file", type=argparse.FileType('r'), required=True)
args = parser.parse_args()

reader = csv.reader(args.input,delimiter="\t")
count_diffs = 0
for line in reader:
    # cols are
    # chrom, pos, ref-allele
    # strain1 allele
    # strain2 allele
    strain1 = line[3]
    strain2 = line[4]
    if strain1 != strain2:
        count_diffs += 1
print(count_diffs)
    