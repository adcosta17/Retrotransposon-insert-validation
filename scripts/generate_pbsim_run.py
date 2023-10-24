import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import random
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Add insertion sequences to contigs')
parser.add_argument('--input', required=True)
parser.add_argument('--output-folder', required=True)
parser.add_argument('--output-prefix', required=True)
parser.add_argument('--output-tsv', required=True)
parser.add_argument('--pbsim-model', required=True)
parser.add_argument('--pbsim-path', required=True)
parser.add_argument('--read-len', required=True, type=int)
args = parser.parse_args()

with open(args.input, 'r') as in_tsv:
    with open(args.output_tsv, 'w') as out_tsv:
        for line in in_tsv:
            row = line.strip().split('\t')
            fa = args.output_folder+"/"+args.output_prefix+"."+str(row[0])+".fa"
            print(args.pbsim_path+"src/pbsim "+fa+" --prefix "+args.output_folder+"/"+args.output_prefix+"."+row[0]+"."+str(args.read_len)+".pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 100 --hmm_model "+args.pbsim_path+args.pbsim_model+" --length-mean "+str(args.read_len)+" --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98 --length-sd "+str(0.2*args.read_len))
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+row[0]+"."+str(args.read_len)+".pbsim_0001.maf")
            print("rm "+args.output_folder+"/"+args.output_prefix+"."+row[0]+"."+str(args.read_len)+".pbsim_0001.ref")
            out_tsv.write(row[0]+"\t"+str(args.read_len)+"\t"+args.output_folder+"/"+args.output_prefix+"."+row[0]+"."+str(args.read_len)+".pbsim_0001.fastq\t"+fa+"\n")

