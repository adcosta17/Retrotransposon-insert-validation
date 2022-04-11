import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import gzip
import random
from intervaltree import Interval, IntervalTree



parser = argparse.ArgumentParser( description='Generate New fastqs with deletions in reads')
parser.add_argument('--mat-ref-inserts', required=True)
parser.add_argument('--pat-ref-inserts', required=True)
parser.add_argument('--mat-contigs', required=True)
parser.add_argument('--pat-contigs', required=True)
parser.add_argument('--out-mat-contigs', required=True)
parser.add_argument('--out-pat-contigs', required=True)
args = parser.parse_args()

# Store positions of hap to ref inserts on the ref. Will then identify reads that align to these positions on the haplotypes end to end
# Will then delete the insert sequence on some fraction of the reads that do align. 
mat_ref_inserts = defaultdict(IntervalTree)
combined_mat_pat_ref_inserts = defaultdict(IntervalTree)
count = 0
with open(args.mat_ref_inserts, 'r') as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        chrom = row[3]
        start = int(row[4])
        end = int(row[5])
        if row[8] != "PASS":
            continue
        mat_ref_inserts[chrom][start:end] = row

#print("pat")
pat_ref_inserts = defaultdict(IntervalTree)
count = 0
with open(args.pat_ref_inserts, 'r') as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        chrom = row[3]
        start = int(row[4])
        end = int(row[5])
        if row[8] != "PASS":
            continue
        pat_ref_inserts[chrom][start:end] = row

seen = {}
# For each contig that has inserts in the mat and pat assembly, remove the inserts
with open(args.out_pat_contigs, 'w') as out_pat, open(args.out_mat_contigs, 'w') as out_mat:
    with pysam.FastaFile(args.pat_contigs) as pat_contigs, pysam.FastaFile(args.mat_contigs) as mat_contigs:
        for chrom in mat_ref_inserts:
            seen[chrom] = 1
            seq = mat_contigs.fetch(chrom)
            pos = 0
            new_seq = ""
            positions = []
            for item in mat_ref_inserts[chrom]:
                start = item.begin
                end = item.end
                positions.append([start, end])
            positions.sort(key=lambda y: y[0])
            for item in positions:
                start = item[0]
                end = item[1]
                new_seq += seq[pos:start]
                pos = end
            new_seq += seq[pos:]
            #print("Old: "+str(len(seq)))
            #print("New: "+str(len(new_seq)))
            out_mat.write(">"+chrom+"\n"+new_seq+"\n")
        for chrom in pat_ref_inserts:
            seen[chrom] = 1
            seq = pat_contigs.fetch(chrom)
            pos = 0
            new_seq = ""
            positions = []
            for item in pat_ref_inserts[chrom]:
                start = item.begin
                end = item.end
                positions.append([start, end])
            positions.sort(key=lambda y: y[0])
            for item in positions:
                start = item[0]
                end = item[1]
                new_seq += seq[pos:start]
                pos = end
            new_seq += seq[pos:]
            #print("Old: "+str(len(seq)))
            #print("New: "+str(len(new_seq)))
            out_pat.write(">"+chrom+"\n"+new_seq+"\n")
    with pysam.FastxFile(args.pat_contigs) as pat_contigs, pysam.FastxFile(args.mat_contigs) as mat_contigs:
        for entry in pat_contigs:
            if entry.name not in seen:
                out_pat.write(">"+entry.name+"\n"+entry.sequence+"\n")
        for entry in mat_contigs:
            if entry.name not in seen:
                out_mat.write(">"+entry.name+"\n"+entry.sequence+"\n")
