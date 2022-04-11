import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import random
from intervaltree import Interval, IntervalTree

random.seed(10)

parser = argparse.ArgumentParser( description='Generate New fastqs with varying coverage levels')
parser.add_argument('--input-fastq', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--coverage', required=True, type=int)
parser.add_argument('--test-number', required=True, type=int)
parser.add_argument('--total-coverage', default=60, type=int)
args = parser.parse_args()
tests = {}
for i in range(args.test_number):
    tests[i+1] = open(args.output_prefix+"_"+str(args.coverage)+"_"+str(i+1)+".fastq", 'w')

count = 0
with pysam.FastxFile(args.input_fastq) as fin:
    for entry in fin:
        count += 1
        if count %1000000 == 0:
            print(count)
        n = random.randint(1, args.total_coverage)
        if n <= (args.coverage*args.test_number):
            # Write to a test sample
            for i in range(args.test_number):
                if n <= ((i+1)*args.coverage):
                    tests[i+1].write("@"+entry.name+"\n")
                    tests[i+1].write(entry.sequence+"\n+\n")
                    tests[i+1].write(entry.quality+"\n")
                    break
        # Otherwise we skip this read if we're downsampling

for i in range(args.test_number):
    tests[i+1].close()