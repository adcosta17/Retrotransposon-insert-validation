import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import random
from intervaltree import Interval, IntervalTree

random.seed(10)

parser = argparse.ArgumentParser( description='Generate New fastqs with spliced in reads')
parser.add_argument('--fastq-folder', default="fastq")
parser.add_argument('--fastq-suffix', default="Guppy_4.0.11_prom.fastq.gz")
parser.add_argument('--input-fastq-1', required=True)
parser.add_argument('--input-fastq-2', required=True)
parser.add_argument('--output-fastq-1', required=True)
parser.add_argument('--output-fastq-2', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--centromeres', required=True)
parser.add_argument('--inserts-tsv', required=True)
args = parser.parse_args()

centromeres = defaultdict(IntervalTree)
with open(args.centromeres, 'r') as in_contromeres:
    count = 0
    for line in in_contromeres:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        chrom = row[1]
        start = int(row[2])
        end = int(row[3])
        key = chrom + ":" + str(start) + "-" + str(end)
        centromeres[chrom][start:end] = key

# Read in combined tsv, store positions of insertions not found in this sample
# Select one of the samples that support an insert at the position, store it and assign a number for which output fastq it is getting sent to
inserts = defaultdict(list)
count = 0
with open(args.inserts_tsv, 'r') as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        pos = row[2]
        chrom = pos.split(':')[0]
        start = int(pos.split(':')[1].split('-')[0])
        end = int(pos.split(':')[1].split('-')[1])
        nearby = centromeres[chrom][start:end]
        if len(nearby) > 0 or row [3] == "Ambiguous":
            continue
        inserts[row[3]].append(row)

# Once we have the reads for each sample open up the files and begin writing
print("Read\tPosition\tInsertLength\tAnnotation\tSourceRead\tSourcePosition\tSeq")
types = ["LINE", "SINE", "SVA", "ERV"]
with open(args.output_fastq_1, 'w') as fout_1, open(args.output_fastq_2, 'w') as fout_2:
    with pysam.FastxFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_3_"+args.fastq_suffix) as fin_1:
        for entry in fin_1:
            n = random.randint(1, 10000)
            if n > 100 or len(entry.sequence) <= 1000:
                fout_1.write("@"+entry.name+"\n")
                fout_1.write(entry.sequence+"\n+\n")
                fout_1.write(entry.quality+"\n")
            else:
                insert_type = types[random.randint(0, 3)]
                row = inserts[insert_type][random.randint(0, len(inserts[insert_type])-1)]
                # Splice in insertion
                insert_seq = row[4]
                insert_pos = random.randint(500, len(entry.sequence)-501)
                seq = entry.sequence[0:insert_pos]+insert_seq+entry.sequence[insert_pos:]
                qual = '='*len(seq)
                fout_1.write("@"+entry.name+"\n")
                fout_1.write(seq+"\n+\n")
                fout_1.write(qual+"\n")
                print(entry.name+"\t"+str(insert_pos)+"\t"+str(len(insert_seq))+"\t"+insert_type+"\t"+row[0]+"\t"+row[2]+"\t"+insert_seq)
    with pysam.FastxFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_5_"+args.fastq_suffix) as fin_2:
        for entry in fin_2:
            n = random.randint(1, 10000)
            if n > 100 or len(entry.sequence) <= 1000:
                fout_2.write("@"+entry.name+"\n")
                fout_2.write(entry.sequence+"\n+\n")
                fout_2.write(entry.quality+"\n")
            else:
                insert_type = types[random.randint(0, 3)]
                row = inserts[insert_type][random.randint(0, len(inserts[insert_type])-1)]
                # Splice in insertion
                insert_seq = row[4]
                insert_pos = random.randint(500, len(entry.sequence)-501)
                seq = entry.sequence[0:insert_pos]+insert_seq+entry.sequence[insert_pos:]
                qual = '='*len(seq)
                fout_2.write("@"+entry.name+"\n")
                fout_2.write(seq+"\n+\n")
                fout_2.write(qual+"\n")
                print(entry.name+"\t"+str(insert_pos)+"\t"+str(len(insert_seq))+"\t"+insert_type+"\t"+row[0]+"\t"+row[2]+"\t"+insert_seq)

