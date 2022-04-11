import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re
import json

def get_match_mismatch(cigarstring):
    count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            count += int(cg[:cg.find("=")])
    return count

def get_read_length(cigarstring):
    # Gets the read length based on the Cigar String
    count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
        elif cg.endswith('='):
            count += int(cg[:cg.find("=")])
    return count

def get_matches(record):
    match_mismatch = get_match_mismatch(record.cigarstring)
    tags = record.get_tags(with_value_type=False)
    mismatches = 0
    for tag in tags:
        if tag[0] == "NM":
            matches = int(tag[1])
            break
    matches = match_mismatch - mismatches
    return matches
    
def get_read_percent_identity(record):
    read_len = get_read_length(record.cigarstring)
    matches = get_matches(record)
    return matches/read_len

def get_aligned_percent_identity(record):
    read_len = record.reference_length
    matches = get_matches(record)
    return matches/read_len

def get_alignment_score(record):
    tags = record.get_tags(with_value_type=False)
    score = 0
    for tag in tags:
        if tag[0] == "AS":
            score = int(tag[1])
            break
    return score

def get_s1_score(record):
    tags = record.get_tags(with_value_type=False)
    score = 0
    for tag in tags:
        if tag[0] == "s1":
            score = int(tag[1])
            break
    return score
        

def get_s1_score(record):
    tags = record.get_tags(with_value_type=False)
    score = 0
    for tag in tags:
        if tag[0] == "s1":
            score = int(tag[1])
            break
    return score


def get_s2_score(record):
    tags = record.get_tags(with_value_type=False)
    score = 0
    for tag in tags:
        if tag[0] == "s2":
            score = int(tag[1])
            break
    return score


def get_de(record):
    tags = record.get_tags(with_value_type=False)
    score = 0
    for tag in tags:
        if tag[0] == "de":
            score = float(tag[1])
            break
    return score

parser = argparse.ArgumentParser( description='Get stats of assemblies mapped to haplotypes for each sample')
parser.add_argument('--sample', required=True)
parser.add_argument('--mapped-to-sample-list', default="HG00438,HG00621,HG00673,HG00735,HG00741,HG01071,HG01106,HG01123")
parser.add_argument('--input-folder', required=True)
parser.add_argument('--min-mapping-qual', type=int, default=60)
parser.add_argument('--min-mapped-fraction', type=float, default=0.3)
parser.add_argument('--min-count', type=int, default=100000)
args = parser.parse_args()

reads_that_map = {}

reads_to_use = {}
with open(args.input_folder+"/"+args.sample+"/"+args.sample+"_chr20_reads.txt", 'r') as in_reads:
    for line in in_reads:
        reads_to_use[line.strip()] = 1

samples_list = args.mapped_to_sample_list.split(',')
for sample in samples_list:
    print(sample)
    # Check each sample
    for hap in ["mat","pat"]:
        print(hap)
        samfile = pysam.AlignmentFile(args.input_folder+"/"+args.sample+"/"+args.sample+".reads_mapped_to."+sample+"."+hap+".bam")
        for record in samfile.fetch():
            if record.query_name not in reads_to_use:
                continue
            if record.query_name not in reads_that_map:
                reads_that_map[record.query_name] = {}
                reads_that_map[record.query_name]["read_length"] = get_read_length(record.cigarstring)
                reads_that_map[record.query_name]["sample_origin"] = args.sample
            if not record.is_secondary and record.mapping_quality > 0 and not record.is_supplementary:
                read_percent_id = get_read_percent_identity(record)
                aligned_percent_id = get_aligned_percent_identity(record)
                score = get_alignment_score(record)
                s1 = get_s1_score(record)
                s2 = get_s2_score(record)
                de = get_de(record)
                reads_that_map[record.query_name][sample+"_"+hap] = {"ref_name": record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end), "alen": record.reference_length, "mapq": record.mapping_quality, "percent_identity": read_percent_id, "aligned_percent_id": aligned_percent_id, "alignment_score": score, "s1_score": s1, "s2_score": s2, "de": de}
print("combined")
samfile = pysam.AlignmentFile(args.input_folder+"/"+args.sample+"/"+args.sample+".reads_mapped_to.combined.bam")
for record in samfile.fetch():
    if record.query_name not in reads_to_use:
        continue
    if record.query_name not in reads_that_map:
        reads_that_map[record.query_name] = {}
        reads_that_map[record.query_name]["read_length"] = get_read_length(record.cigarstring)
        reads_that_map[record.query_name]["sample_origin"] = args.sample
    if not record.is_secondary and record.mapping_quality > 0 and not record.is_supplementary:
        read_percent_id = get_read_percent_identity(record)
        aligned_percent_id = get_aligned_percent_identity(record)
        score = get_alignment_score(record)
        s1 = get_s1_score(record)
        s2 = get_s2_score(record)
        de = get_de(record)
        reads_that_map[record.query_name]["combined"] = {"ref_name": record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end), "alen": record.reference_length, "mapq": record.mapping_quality, "percent_identity": read_percent_id, "aligned_percent_id": aligned_percent_id, "alignment_score": score, "s1_score": s1, "s2_score": s2, "de": de}
print("grch38")
samfile = pysam.AlignmentFile(args.input_folder+"/"+args.sample+"/"+args.sample+".reads_mapped_to.grch38.bam")
for record in samfile.fetch():
    if record.query_name not in reads_to_use:
        continue
    if record.query_name not in reads_that_map:
        reads_that_map[record.query_name] = {}
        reads_that_map[record.query_name]["read_length"] = get_read_length(record.cigarstring)
        reads_that_map[record.query_name]["sample_origin"] = args.sample
    if not record.is_secondary and record.mapping_quality > 0 and not record.is_supplementary:
        read_percent_id = get_read_percent_identity(record)
        aligned_percent_id = get_aligned_percent_identity(record)
        score = get_alignment_score(record)
        s1 = get_s1_score(record)
        s2 = get_s2_score(record)
        de = get_de(record)
        reads_that_map[record.query_name]["grch38"] = {"ref_name": record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end), "alen": record.reference_length, "mapq": record.mapping_quality, "percent_identity": read_percent_id, "aligned_percent_id": aligned_percent_id, "alignment_score": score, "s1_score": s1, "s2_score": s2, "de": de}
                

with open(args.sample+".json", 'w') as fp:
    json.dump(reads_that_map, fp)

# Check combined and grch38
