import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re

def get_read_len(cigarstring):
    # Count up the position on the read until we get to a deletion
    count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        if cg.endswith('='):
            count += int(cg[:cg.find("=")])
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
    return count

parser = argparse.ArgumentParser( description='Filter Reads based on their alignment to each assembly')
parser.add_argument('--mat-bam', required=True)
parser.add_argument('--pat-bam', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--input-folder', required=True)
parser.add_argument('--min-mapping-qual', type=int, default=60)
parser.add_argument('--min-mapped-fraction', type=float, default=0.3)
args = parser.parse_args()

mat_contigs = defaultdict(int)
pat_contigs = defaultdict(int)
with pysam.FastxFile(args.input_folder+'/'+args.sample+"/"+args.sample+"_chr20.mat.assembly.fa") as in_mat:
	for entry in in_mat:
		mat_contigs[entry.name.split('_')[1]] = len(entry.sequence)
with pysam.FastxFile(args.input_folder+'/'+args.sample+"/"+args.sample+"_chr20.pat.assembly.fa") as in_pat:
	for entry in in_pat:
		pat_contigs[entry.name.split('_')[1]] = len(entry.sequence)

reads_to_use = {}
sam_reader = pysam.AlignmentFile(args.mat_bam)
for contig in mat_contigs:
	try:
		mappings = sam_reader.fetch(contig)
	except:
		continue
	for record in mappings:
		r_len = get_read_len(record.cigarstring)
		if record.mapping_quality >= args.min_mapping_qual and (record.query_alignment_end - record.query_alignment_start)/r_len > args.min_mapped_fraction:
			reads_to_use[record.query_name] = 1

sam_reader = pysam.AlignmentFile(args.pat_bam)
for contig in pat_contigs:
	try:
		mappings = sam_reader.fetch(contig)
	except:
		continue
	for record in mappings:
		r_len = get_read_len(record.cigarstring)
		if record.mapping_quality >= args.min_mapping_qual and (record.query_alignment_end - record.query_alignment_start)/r_len > args.min_mapped_fraction:
			reads_to_use[record.query_name] = 1

for read in reads_to_use:
	print(read)
