import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re
import json
from intervaltree import Interval, IntervalTree
from statistics import mean

def match(name, truth_reads):
    for read in truth_reads:
        if read in name:
            return read
    return None

def same_family(family, annotation):
    if "line" in family.lower() and "line" in annotation.lower():
        return True
    elif "alu" in family.lower() and "alu" in annotation.lower():
        return True
    elif "sva" in family.lower() and "sva" in annotation.lower():
        return True
    else:
        return False

def get_nearby(chrom, start, end, seen):
    if chrom in seen:
        count = 0
        i = start
        while i < end:
            if i in seen[chrom]:
                count += seen[chrom][i]
            i += 1
        if count >= 1:
            return True
        return False
    return True

def get_nearby_tp(chrom, start, end, seen):
    if chrom in seen:
        count = 0
        i = start
        while i < end:
            if i in seen[chrom]:
                count += seen[chrom][i]
            i += 1
        if count > 0:
            return True
        return False
    return True



def get_pipeline_stats(file, window, truth_reads, truth_dict, min_size, max_size):
    base_tree = defaultdict(IntervalTree)
    fp_dict = defaultdict(IntervalTree)
    seen = {}
    all_inserts = {}
    tp = 0
    fn = 0
    fp = 0
    found_count = defaultdict(int)
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
        all_inserts[chrom] = defaultdict(int)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            if row[3] in truth_reads:
                found_count[int(truth_reads[row[3]][0])] = 1
            if "PASS" not in line:
                continue
            if "ambiguous" in line:
                continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:end+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(item.data[0], row[9]):
                        seen[chrom][int(item.begin)] = 1
            if chrom not in all_inserts:
                all_inserts[chrom] = defaultdict(int)
            all_inserts[chrom][start] = 1
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line:
                continue
            if "ambiguous" in line:
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:end+window]
            found = False
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(item.data[0], row[9]):
                        found = True
            if found:
                continue
            nearby = get_nearby(chrom, start-50, end+50, all_inserts)
            nearby_fp = fp_dict[chrom][start:end]
            if nearby and len(nearby_fp) == 0:
                fp_dict[chrom][(start-50):(end+50)] = start
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.begin) not in seen[chrom]:
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
    for chrom in all_inserts:
        for pos in fp_dict[chrom]:
            fp += 1
    missing_count = 0
    for i in range(500):
        if found_count[i] == 0:
            missing_count += 1
    return tp, fn, fp, missing_count


def get_somrit_stats(file, window, truth_reads, truth_dict, min_size, max_size):
    base_tree = defaultdict(IntervalTree)
    fp_dict = defaultdict(IntervalTree)
    seen = {}
    all_inserts = {}
    tp = 0
    fn = 0
    fp = 0
    found_count = defaultdict(int)
    for chrom in truth_dict:
        seen[chrom] = defaultdict(int)
        all_inserts[chrom] = defaultdict(int)
    count = 0
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            row = line.strip().split('\t')
            for item in row[6].split(','):
                if item.split(':')[0] == "Sample3" and item.split(':')[1] in truth_reads:
                    found_count[int(truth_reads[item.split(':')[1]][0])] += 1
            if "PASS" not in line:
                continue
            if "PossiblyNovel_NotInAssembly" not in line:
                continue
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:end+window]
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(item.data[0], row[8]):
                        seen[chrom][int(item.begin)] = 1
            if chrom not in all_inserts:
                all_inserts[chrom] = defaultdict(int)
            all_inserts[chrom][start] = 1
    with open(file) as in_calls:
        for line in in_calls:
            if count == 0:
                count += 1
                continue
            if "PASS" not in line:
                continue
            if "PossiblyNovel_NotInAssembly" not in line:
                continue
            row = line.strip().split('\t')
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            # Check if the position is within window size of a truth call
            nearby = truth_dict[chrom][start-window:end+window]
            found = False
            if len(nearby) > 0:
                # Have a tp
                for item in nearby:
                    #print(item.data)
                    if same_family(item.data[0], row[8]):
                        found = True
            if found:
                continue
            nearby = get_nearby(chrom, start-50, end+50, all_inserts)
            nearby_fp = fp_dict[chrom][start:end]
            if nearby and len(nearby_fp) == 0:
                fp_dict[chrom][(start-50):(end+50)] = start
    for chrom in truth_dict:
        for item in truth_dict[chrom]:
            if int(item.begin) not in seen[chrom]:
                fn += 1
    for chrom in seen:
        for pos in seen[chrom]:
            tp += 1
    for chrom in all_inserts:
        for pos in fp_dict[chrom]:
            fp += 1
    missing_count = 0
    for i in range(500):
        if found_count[i] == 0:
            missing_count += 1
    return tp, fn, fp, missing_count




parser = argparse.ArgumentParser( description='Get metrics for xtea and tldr before and after runs')
parser.add_argument('--inserts-tsv', required=True)
parser.add_argument('--spiked-reads', required=True)
parser.add_argument('--pipeline', required=True)
parser.add_argument('--somrit', required=True)
parser.add_argument('--min-size', default=0, type=int)
parser.add_argument('--max-size', default=10000, type=int)
parser.add_argument('--window-size', type=int, default=1000)
parser.add_argument('--chrom-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
args = parser.parse_args()


# Start by reading in the list of spiked in reads and the positions they come from
# Then read in the metadata on those positions that will tell us more about the positions
spike_in_reads = defaultdict(list)
spike_in_positions = defaultdict(IntervalTree)
insert_data = defaultdict(list)
with open(args.inserts_tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        insert_data[row[0]] = row
        spike_in_positions[row[5]][int(row[6]):int(row[6])+1] = [row[7], len(row[8])]

# Get the per read stats for somrit and tldr
# Positions for sniffles and xt
with open(args.spiked_reads, 'r') as in_reads:
    for line in in_reads:
        row = line.strip().split('\t')
        values = [row[0], insert_data[row[0]][5], int(insert_data[row[0]][6]), insert_data[row[0]][7], insert_data[row[0]][0], insert_data[row[0]][1], insert_data[row[0]][2], insert_data[row[0]][3], insert_data[row[0]][4]]
        if len(insert_data[row[0]][8]) < args.min_size or len(insert_data[row[0]][8]) > args.max_size:
            continue
        spike_in_reads[row[1]] = values

 # Get stats for each case
tp_somrit, fn_somrit, fp_somrit, found_somrit_count  = get_somrit_stats(args.somrit,args.window_size,spike_in_reads,spike_in_positions, args.min_size, args.max_size)
tp_pipeline, fn_pipeline, fp_pipeline, found_pipeline_count = get_pipeline_stats(args.pipeline,args.window_size,spike_in_reads,spike_in_positions, args.min_size, args.max_size)
som_prec = 0
if tp_somrit + fp_somrit > 0:
    som_prec = float(tp_somrit)/(tp_somrit+fp_somrit)
som_rec = 0
if tp_somrit + fn_somrit > 0:
    som_rec = float(tp_somrit)/(tp_somrit+fn_somrit)
pipe_rec = 0
if tp_pipeline + fn_pipeline > 0:
    pipe_rec = float(tp_pipeline)/(tp_pipeline+fn_pipeline)
pipe_prec = 0
if tp_pipeline + fp_pipeline > 0:
    pipe_prec = float(tp_pipeline)/(tp_pipeline+fp_pipeline)

to_print = [str(tp_somrit), str(fn_somrit), str(fp_somrit), str(som_rec), str(som_prec), str(tp_pipeline), str(fn_pipeline), str(fp_pipeline), str(pipe_rec), str(pipe_prec), str(found_somrit_count), str(found_pipeline_count)]
print("\t".join(to_print))

