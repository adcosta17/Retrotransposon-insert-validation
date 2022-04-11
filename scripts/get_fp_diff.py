import pysam
import argparse
import sys
import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import gzip
import re
import math

parser = argparse.ArgumentParser( description='Filter Reads based on their alignment to each assembly')
parser.add_argument('--before', required=True)
parser.add_argument('--after', required=True)
args = parser.parse_args()

before = {}
with open(args.before) as in_before:
    for line in in_before:
        row = line.strip().split("\t")
        start = str(math.ceil(int(row[1])/1000)*1000)
        end = str(math.ceil(int(row[2])/1000)*1000)
        #before[row[0]+":"+start+"-"+end] = row
        before[row[0]+":"+row[3]] = row

after = {}
with open(args.after) as in_after:
    for line in in_after:
        row = line.strip().split("\t")
        start = str(math.ceil(int(row[1])/1000)*1000)
        end = str(math.ceil(int(row[2])/1000)*1000)
        #after[row[0]+":"+start+"-"+end] = row
        after[row[0]+":"+row[3]] = row

in_both = {}
only_before = {}
only_after = {}
for pos in after:
    if pos in before:
        in_both[pos] = 1
        print("Both\t"+"\t".join(after[pos]))
    else:
        only_after[pos] = 1
        print("After\t"+"\t".join(after[pos]))

for pos in before:
    if pos not in after:
        only_before[pos] = 1
        print("Before\t"+"\t".join(before[pos]))



