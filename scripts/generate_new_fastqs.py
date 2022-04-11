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
parser.add_argument('--suffix', default=".insert_supporting_reads.tsv")
parser.add_argument('--folder', default="reads_assembly_mapped")
parser.add_argument('--fastq-folder', default="fastq")
parser.add_argument('--fastq-suffix', default="Guppy_4.0.11_prom.fastq.gz")
parser.add_argument('--input-fastq-1', required=True)
parser.add_argument('--input-fastq-2', required=True)
parser.add_argument('--output-fastq-1', required=True)
parser.add_argument('--output-fastq-2', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--centromeres', required=True)
parser.add_argument('--all-samples', required=True)
parser.add_argument('--combined-tsv', required=True)
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
positions = {}
count = 0
with open(args.combined_tsv, 'r') as in_combined:
    for line in in_combined:
        if args.sample in line:
            continue
        row = line.strip().split('\t')
        chrom = row[0].split(':')[0]
        start = int(row[0].split(':')[1].split('-')[0])
        end = int(row[0].split(':')[1].split('-')[1])
        nearby = centromeres[chrom][start:end]
        if len(nearby) > 0:
            continue
        # Only consider positions that are PASS. Need the insert to be annotated to a repeat
        if row[2] == "FAIL":
            continue
        sample = random.choice(row[1].split(',')).split('_')[0]
        if sample not in positions:
            positions[sample] = {}
        positions[sample][row[0]] = count % 2
        count += 1

# For each sample in positions, Read in its position file, and get the read names at each position for that sample
one_reads = {}
two_reads = {}
for sample in args.all_samples.split(','):
    if sample == args.sample:
        continue
    one_reads[sample] = defaultdict(list)
    two_reads[sample] = defaultdict(list)
    with open(sample+"/"+args.folder+"/"+sample+args.suffix, 'r') as in_reads:
        count = 0
        for line in in_reads:
            if count == 0:
                count = 1
                continue
            row = line.strip().split('\t')
            if row[2] in positions[sample]:
                if positions[sample][row[2]] == 0:
                    one_reads[sample][row[2]].append(row)
                else:
                    two_reads[sample][row[2]].append(row)

# Once we have the reads for each sample open up the files and begin writing
print("Read\tSource\tPosition\tInsertLength\tAnnotation")
s3_total_len = 0
s5_total_len = 0
s3_cov = 0
s5_cov = 0
with open(args.output_fastq_1, 'w') as fout_1, open(args.output_fastq_2, 'w') as fout_2:
    with pysam.FastxFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_3_"+args.fastq_suffix) as fin_1:
        for entry in fin_1:
            fout_1.write("@"+entry.name+"\n")
            fout_1.write(entry.sequence+"\n+\n")
            fout_1.write(entry.quality+"\n")
            s3_total_len += len(entry.sequence)
    with pysam.FastxFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_5_"+args.fastq_suffix) as fin_2:
        for entry in fin_2:
            fout_2.write("@"+entry.name+"\n")
            fout_2.write(entry.sequence+"\n+\n")
            fout_2.write(entry.quality+"\n")
            s5_total_len += len(entry.sequence)
    s3_cov = int(s3_total_len/3200000000)
    s5_cov = int(s5_total_len/3200000000)
    for sample in args.all_samples.split(','):
        if sample == args.sample:
            continue
        with pysam.FastaFile(sample+"/"+args.fastq_folder+"/"+sample+"_1_"+args.fastq_suffix) as fin:
            seen = {}
            for position in one_reads[sample]:
                n = random.randint(1, min(len(one_reads[sample][position]), s3_cov))
                rows = random.sample(one_reads[sample][position], n)
                for row in rows:
                    if row[0] not in seen:
                        read_to_use = row[0]
                        seq = fin.fetch(read_to_use)
                        fout_1.write("@"+read_to_use+"\n"+seq+"\n+\n"+'='*len(seq)+'\n')
                        print(read_to_use+"\t"+sample+"\t"+position+"\t"+row[1]+"\t"+row[3])
                        seen[read_to_use] = 1
            for position in two_reads[sample]:
                n = random.randint(1, min(len(two_reads[sample][position]), s5_cov))
                rows = random.sample(two_reads[sample][position], n)
                for row in rows:
                    if row[0] not in seen:
                        read_to_use = row[0]
                        seq = fin.fetch(read_to_use)
                        fout_2.write("@"+read_to_use+"\n"+seq+"\n+\n"+'='*len(seq)+'\n')
                        print(read_to_use+"\t"+sample+"\t"+position+"\t"+row[1]+"\t"+row[3])
                        seen[read_to_use] = 1

