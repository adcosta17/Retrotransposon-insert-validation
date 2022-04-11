import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re

def check_no_indels(cg_tag):
    cigar = cg_tag.split(':')[2]
    for cg in re.findall('[0-9]*[A-Z]', cigar):
        if cg.endswith('I'):
            if int(cg[:cg.find("I")]) >= 100:
                return False
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) >= 100:
                return False
    return True

def get_positions(cg_tag, asm_start, asm_end, read_to_contig_start):
    # Gets the insertion start and end  based on the Cigar String
    cigar = cg_tag.split(':')[2]
    read_count = 0
    ref_count = read_to_contig_start
    read_start = 0
    read_end = 0
    #print(str(asm_start)+ " "+ str(asm_end))
    for cg in re.findall('[0-9]*[A-Z=]', cigar):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            ref_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('P'):
            read_count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        if ref_count >= asm_start and read_start == 0:
            if cg.endswith('M'):
                read_start = read_count - int(cg[:cg.find("M")])
            elif cg.endswith('X'):
                read_start = read_count - int(cg[:cg.find("X")])
            else:
                read_start = read_count
        if ref_count >= asm_end and read_end == 0:
                read_end = read_count
    #print(read_count)
    return read_start, read_end

def get_family(annotation):
    if "LINE" in annotation:
        return "LINE"
    if "SINE" in annotation:
        return "SINE"
    if "ERV" in annotation:
        return "ERV"
    if "SVA" in annotation:
        return "SVA"
    return "Ambiguous"

def check_pos_for_line(line, min_mapping_qual, passing_inserts, insert_annoations):
    row = line.decode().strip().split('\t')
    if int(row[11]) <= min_mapping_qual:
        return
    if row[5] in passing_inserts:
        for pos in passing_inserts[row[5]]:
            asm_start = int(pos.split('-')[0])
            asm_end = int(pos.split('-')[1])
            read_to_contig_start = int(row[7])
            read_to_contig_end = int(row[8])
            if asm_start - read_to_contig_start >= 500 and read_to_contig_end - asm_end >= 500:
                # Have an overlapping alignmnet. Check the CIGAR string to make sure that there are no large indels
                if check_no_indels(row[22]):
                    # Now have to extract the insertions from the read, parse the CIGAR of the alignment between the read and ref to get the read start and end of the insert
                    read_start, read_end = get_positions(row[22], asm_start, asm_end, read_to_contig_start)
                    try:
                        seq = fq.fetch(row[0])
                        if row[4] == '-':
                            # Reverse complement seq
                            tmp = len(seq) - read_start
                            read_start = len(seq) - read_end
                            read_end = tmp
                        #print(read_start)
                        #print(read_end)
                        if read_start == 0:
                            continue
                        print(row[0]+'\t'+str(asm_end - asm_start)+'\t'+passing_inserts[row[5]][pos]+'\t'+insert_annoations[row[5]][pos]+"\t"+row[5]+":"+pos+"\t"+seq[read_start:read_end+1])
                    except KeyError:
                        #print("Error for "+str(row))
                        pass

parser = argparse.ArgumentParser( description='Parse read alignment to assembly and flag those that align end to end')
parser.add_argument('--input-tsv', required=True)
parser.add_argument('--input-mat-paf', required=True)
parser.add_argument('--input-pat-paf', required=True)
parser.add_argument('--input-fastq', required=True)
parser.add_argument('--min-mapping-qual', type=int, default=50)
args = parser.parse_args()

# First read in tsv of passing insertions
# Get contig and positions on contig for each passing insert
passing_inserts = {}
insert_annoations = {}
with open(args.input_tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0:
            count = 1
            continue
        if "PASS" in row[8]:
            if row[3] not in passing_inserts:
                passing_inserts[row[3]] = defaultdict(str)
            if row[3] not in insert_annoations:
                insert_annoations[row[3]] = defaultdict(str)
            passing_inserts[row[3]][row[4]+"-"+row[5]] = row[0]+":"+row[1]+"-"+row[2]
            insert_annoations[row[3]][row[4]+"-"+row[5]] = get_family(row[9])

# Open paf.gz file
# Read in alignments and only look at reads that map with a mapq60
# For any read that maps through a passing asssembly insertion position, 
# output the read name, position on read where we expect the insertion to be, the type of insert
# also output the genomic position of where we expect it to be based off off the assembly alignment
# We will refer back to this list when generating the spliced in samples
print("Read\tExpectedInsertSize\tLocation\tAnnotation\tAssemblyPosition\tInsertionSequence")
with gzip.open(args.input_mat_paf,'r') as in_mat_paf:
    with gzip.open(args.input_pat_paf,'r') as in_pat_paf:
        with pysam.FastaFile(filename=args.input_fastq) as fq:
            for line in in_mat_paf:
                check_pos_for_line(line, args.min_mapping_qual, passing_inserts, insert_annoations)
            for line in in_pat_paf:
                check_pos_for_line(line, args.min_mapping_qual, passing_inserts, insert_annoations)

