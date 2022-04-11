import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases


def get_id(header, name):
    i = 0
    while i < len(header['SQ']):
        if header['SQ'][i]['SN'] == name:
            return i
        else:
            i += 1
    return 0


parser = argparse.ArgumentParser( description='Convert Paf to Bam')
parser.add_argument('--input-paf', required=True)
parser.add_argument('--input-fastq', required=True)
parser.add_argument('--output-bam', required=True)
args = parser.parse_args()

header = { 'HD': {'VN': '1.0'}, 'SQ': [] }

contigs = defaultdict(int)
with gzip.open(args.input_paf,'r') as in_paf:
    for line in in_paf:
        row = line.decode().strip().split('\t')
        if row[5] not in contigs:
            contigs[row[5]] = int(row[6])
            header['SQ'].append({'LN': int(row[6]), 'SN': row[5]})

# Open paf.gz file
read_count = defaultdict(int)
with gzip.open(args.input_paf,'r') as in_paf:
    with pysam.FastaFile(filename=args.input_fastq) as fq:
        with pysam.AlignmentFile(args.output_bam, "wb", header=header) as outf:
            for line in in_paf:
                row = line.decode().strip().split('\t')
                a = pysam.AlignedSegment()
                a.query_name = row[0]
                if read_count[row[0]] == 0:
                    a.flag = 0
                else:
                    a.flag = 256
                if row[4] == '-':
                    a.flag += 16
                a.reference_id = get_id(header, row[5])
                a.reference_start = int(row[7])
                a.mapping_quality = int(row[11])
                cigar = ""
                for item in row:
                    if item.startswith("cg:Z"):
                        cigar = item
                        break
                a.cigarstring = cigar.split(':')[2]
                seq = fq.fetch(row[0])[int(row[2]):int(row[3])]
                if row[4] == '-':
                    seq = reverse_complement(seq)
                a.query_sequence = seq
                outf.write(a)

