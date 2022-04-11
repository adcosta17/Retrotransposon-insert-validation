import pysam
import argparse
import sys
import csv
from collections import defaultdict
import gzip
import re
import random
from intervaltree import Interval, IntervalTree

def get_read_length(cigarstring):
    contig_count = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            contig_count += int(cg[:cg.find("M")])
        elif cg.endswith('='):
            contig_count += int(cg[:cg.find("=")])
        elif cg.endswith('X'):
            contig_count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            contig_count += int(cg[:cg.find("I")])
        elif cg.endswith('S'):
            contig_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            contig_count += int(cg[:cg.find("H")])
    return contig_count

def get_contig_position(cigarstring, ref_start, pos, is_reverse, q_length):
    contig_count = 0
    #print(pos)
    #print(ref_start)
    ref_count = ref_start
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            contig_count += int(cg[:cg.find("M")])
            ref_count += int(cg[:cg.find("M")])
        elif cg.endswith('='):
            contig_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('X'):
            contig_count += int(cg[:cg.find("X")])
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            contig_count += int(cg[:cg.find("I")])
        elif cg.endswith('S'):
            contig_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            contig_count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('N'):
            ref_count += int(cg[:cg.find("N")])
        if ref_count >= pos:
            if is_reverse:
                return q_length - contig_count
            else:
                return contig_count
    if is_reverse:
        return q_length - contig_count
    else:
        return contig_count

random.seed(100)

parser = argparse.ArgumentParser( description='Generate New fastqs with spliced in reads')
parser.add_argument('--fastq-folder', default="fastq")
parser.add_argument('--fastq-suffix', default="Guppy_4.0.11_prom.fastq.gz")
parser.add_argument('--input-fastq-1', required=True)
parser.add_argument('--input-fastq-2', required=True)
parser.add_argument('--mat-input-bam-1', required=True)
parser.add_argument('--mat-input-bam-2', required=True)
parser.add_argument('--mat-assembly-bam', required=True)
parser.add_argument('--mat-assembly-fa', required=True)
parser.add_argument('--pat-input-bam-1', required=True)
parser.add_argument('--pat-input-bam-2', required=True)
parser.add_argument('--pat-assembly-bam', required=True)
parser.add_argument('--pat-assembly-fa', required=True)
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

contigs_per_sample = {}
contigs_to_use = {}
contigs_per_sample["mat"] = defaultdict(list)
contigs_to_use["mat"] = {}
contigs_per_sample["pat"] = defaultdict(list)
contigs_to_use["pat"] = {}

contig_positions_to_avoid = defaultdict(IntervalTree)
contig_mapped_positions = defaultdict(IntervalTree)

# Read in Assembly to ref bam, Get list of contigs that have a mapq 60 alignment
mat_sam_reader = pysam.AlignmentFile(args.mat_assembly_bam)
for record in mat_sam_reader.fetch():
    contigs_per_sample["mat"][record.query_name].append((record.query_alignment_start, record.query_alignment_end-record.query_alignment_start))
    if record.mapping_quality < 60:
        continue
    q_start = record.query_alignment_start
    q_end = record.query_alignment_end
    if q_start < q_end:
        contig_mapped_positions[record.query_name][q_start:q_end]
        #print(record.query_name+" mapped from "+str(q_start)+" "+str(q_end))
    else:
        contig_mapped_positions[record.query_name][q_end:q_start]
        #print(record.query_name+" mapped from "+str(q_start)+" "+str(q_end))
    # Check centromeric alignments. If part of the contig is aligned to the centromere, mark the positions
    nearby = centromeres[chrom][record.reference_start:record.reference_end]
    if len(nearby) > 0:
        # Get the position on the contig that overlap with the centromere
        for item in nearby:
            start = item.begin
            end = item.end
            q_start = record.query_alignment_start
            q_end = record.query_alignment_end
            q_length = get_read_length(record.cigarstring)
            if record.is_reverse:
                #print("reverse")
                tmp = q_start
                q_start = q_end
                q_end = tmp
            if start > record.reference_start and end < record.reference_end:
                # Centromere is in middle
                #print("1\t"+str(record.reference_start)+"\t"+str(start)+"\t"+str(end), flush=True)
                contig_start = get_contig_position(record.cigarstring, record.reference_start, start, record.is_reverse, q_length)
                contig_end  = get_contig_position(record.cigarstring, record.reference_start, end, record.is_reverse, q_length)
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
            elif start < record.reference_start and end < record.reference_end:
                # overlap at the start of the contig
                #print("2\t"+str(record.reference_start)+"\t"+str(end), flush=True)
                contig_start = q_start
                contig_end  = get_contig_position(record.cigarstring, record.reference_start, end, record.is_reverse, q_length)
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
            elif end > record.reference_end and start > record.reference_start:
                # overlap at end of contig
                #print("3\t"+str(record.reference_start)+"\t"+str(start), flush=True)
                contig_start = get_contig_position(record.cigarstring, record.reference_start, start, record.is_reverse, q_length)
                contig_end = q_end
                #print("3\t"+str(q_start)+"\t"+str(q_end)+"\t"+str(contig_start)+"\t"+str(contig_end), flush=True)
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
            elif start < record.reference_start and end > record.reference_end:
                contig_start = q_start
                contig_end = q_end
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
    contigs_to_use["mat"][record.query_name] = 1
pat_sam_reader = pysam.AlignmentFile(args.pat_assembly_bam)
for record in pat_sam_reader.fetch():
    contigs_per_sample["pat"][record.query_name].append((record.query_alignment_start, record.query_alignment_end-record.query_alignment_start))
    if record.mapping_quality < 60:
        continue
    q_start = record.query_alignment_start
    q_end = record.query_alignment_end
    if q_start < q_end:
        contig_mapped_positions[record.query_name][q_start:q_end]
        #print(record.query_name+" mapped from "+str(q_start)+" "+str(q_end))
    else:
        contig_mapped_positions[record.query_name][q_end:q_start]
        #print(record.query_name+" mapped from "+str(q_start)+" "+str(q_end))
    # Check centromeric alignments. If part of the contig is aligned to the centromere, mark the positions
    nearby = centromeres[chrom][record.reference_start:record.reference_end]
    if len(nearby) > 0:
        # Get the position on the contig that overlap with the centromere
        for item in nearby:
            start = item.begin
            end = item.end
            q_start = record.query_alignment_start
            q_end = record.query_alignment_end
            q_length = get_read_length(record.cigarstring)
            if record.is_reverse:
                tmp = q_start
                q_start = q_end
                q_end = tmp
            if start > record.reference_start and end < record.reference_end:
                # Centromere is in middle
                #print("1\t"+str(record.reference_start)+"\t"+str(start)+"\t"+str(end), flush=True)
                contig_start = get_contig_position(record.cigarstring, record.reference_start, start, record.is_reverse, q_length)
                contig_end  = get_contig_position(record.cigarstring, record.reference_start, end, record.is_reverse, q_length)
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
            elif start < record.reference_start and end < record.reference_end:
                # overlap at the start of the contig
                #print("2\t"+str(record.reference_start)+"\t"+str(end), flush=True)
                contig_start = q_start
                contig_end  = get_contig_position(record.cigarstring, record.reference_start, end, record.is_reverse, q_length)
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
            elif end > record.reference_end and start > record.reference_start:
                # overlap at end of contig
                #print("3\t"+str(record.reference_start)+"\t"+str(start), flush=True)
                contig_start = get_contig_position(record.cigarstring, record.reference_start, start, record.is_reverse, q_length)
                contig_end = q_end
                #print("3\t"+str(contig_start)+"\t"+str(contig_end), flush=True)
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
            elif start < record.reference_start and end > record.reference_end:
                contig_start = q_start
                contig_end = q_end
                if contig_start < contig_end:
                    contig_positions_to_avoid[record.query_name][contig_start:contig_end] = 1
                else:
                    contig_positions_to_avoid[record.query_name][contig_end:contig_start] = 1
    contigs_to_use["pat"][record.query_name] = 1

contig_lengths = defaultdict(int)
pat_contigs = {}
mat_contigs = {}
with pysam.FastxFile(args.mat_assembly_fa) as mat_asm_fa:
        for entry in mat_asm_fa:
            contig_lengths[entry.name] = len(entry.sequence)
            mat_contigs[entry.name] = 1
with pysam.FastxFile(args.pat_assembly_fa) as pat_asm_fa:
        for entry in pat_asm_fa:
            contig_lengths[entry.name] = len(entry.sequence)
            pat_contigs[entry.name] = 1

contigs_to_consider = []
for contig in contigs_to_use["mat"]:
    positions = sorted(contigs_per_sample["mat"][contig])
    count = 0
    pos = 0
    for item in positions:
        start = item[0]
        length = item[1]
        if start+length <= pos:
            # This position has already been aligned
            continue
        if start <= pos:
            # have an overlapping region 
            count = count+(start+length-pos)
            pos = start+length
        else:
            # have a gap
            count = count + length
            pos = start+length
    if contig_lengths[contig] > 0:
        if count/contig_lengths[contig] > 0.3:
            contigs_to_consider.append(contig)

for contig in contigs_to_use["pat"]:
    positions = sorted(contigs_per_sample["pat"][contig])
    count = 0
    pos = 0
    for item in positions:
        start = item[0]
        length = item[1]
        if start+length <= pos:
            # This position has already been aligned
            continue
        if start <= pos:
            # have an overlapping region 
            count = count+(start+length-pos)
            pos = start+length
        else:
            # have a gap
            count = count + length
            pos = start+length
    if contig_lengths[contig] > 0:
        if count/contig_lengths[contig] > 0.3:
            contigs_to_consider.append(contig)

# Read in combined tsv, store positions of insertions not found in this sample
# Select one of the samples that support an insert at the position, store it and assign a number for which output fastq it is getting sent to
inserts = defaultdict(list)
insert_contig_pos = defaultdict(IntervalTree)
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
        contig = row[4].split(':')[0]
        c_start = int(row[4].split(':')[1].split('-')[0])
        c_end = int(row[4].split(':')[1].split('-')[1])
        if len(row) > 5:
            inserts[row[4]].append(row)
            insert_contig_pos[contig][c_start:c_end] = 1

def get_insert_read_pos(cigar, ref_start, asm_pos):
    # Gets the insertion start and end  based on the Cigar String
    read_count = 0
    ref_count = ref_start
    read_end = 0
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
        if ref_count >= asm_pos:
            if cg.endswith('M') or cg.endswith('X') or  cg.endswith('='):
                return read_count - (asm_pos - ref_count)
            else:
                return read_count
    return -1

# Randomly select positions from the alignment of reads to the assembly that are distant from the insertions
count = 0
mat_sam_reader_1 = pysam.AlignmentFile(args.mat_input_bam_1)
mat_sam_reader_2 = pysam.AlignmentFile(args.mat_input_bam_2)
pat_sam_reader_1 = pysam.AlignmentFile(args.pat_input_bam_1)
pat_sam_reader_2 = pysam.AlignmentFile(args.pat_input_bam_2)
reads = {}
reads[1] = defaultdict(list)
reads[2] = defaultdict(list)
while count < 2000:
    # Insert 2000 positions
    n = random.randint(0, 1)
    if n == 0:
        # input 1
        contig = random.choice(contigs_to_consider)
        pos = random.randint(1000, contig_lengths[contig] - 1000)
        if len(insert_contig_pos[contig][pos-1000:pos+1000]) > 0:
            continue
        if len(contig_positions_to_avoid[contig][pos:pos+1]) > 0:
            continue
        #if len(contig_mapped_positions[contig][pos:pos+1]) == 0:
        #    continue
        #print("Found: "+contig+":"+str(pos))
        # Position passes, now have to get reads
        try:
            read_count = 0
            if contig in pat_contigs:
                sam_reader_1 = pat_sam_reader_1
            else:
                sam_reader_1 = mat_sam_reader_1
            for record in sam_reader_1.fetch(region=contig+":"+str(pos)+"-"+str(pos+1)):
                # Get the read positions of where pos is
                aligned_frac = (record.query_alignment_end - record.query_alignment_start)/get_read_length(record.cigarstring)
                if record.mapq < 60 or aligned_frac < 0.7:
                    continue
                read_count += 1
                read_pos = get_insert_read_pos(record.cigarstring, record.reference_start, pos)
                reads[1][contig+":"+str(pos)].append([record.query_name, read_pos])
            if read_count < 5:
                continue
            insert_contig_pos[contig][pos:pos+1] = 1
            count += 1
        except ValueError:
            pass

    else:
        # input 2
        contig = random.choice(contigs_to_consider)
        pos = random.randint(1000, contig_lengths[contig] - 1000)
        if len(insert_contig_pos[contig][pos-1000:pos+1000]) > 0:
            continue
        if len(contig_positions_to_avoid[contig][pos:pos+1]) > 0:
            continue
        #if len(contig_mapped_positions[contig][pos:pos+1]) == 0:
            #print(contig+" Not mapped at "+str(pos))
        #    continue
        #print("Found: "+contig+":"+str(pos))
        # Position passes, now have to get reads
        try:
            read_count = 0
            if contig in pat_contigs:
                sam_reader_2 = pat_sam_reader_2
            else:
                sam_reader_2 = mat_sam_reader_2
            for record in sam_reader_2.fetch(region=contig+":"+str(pos)+"-"+str(pos+1)):
                # Get the read positions of where pos is
                aligned_frac = (record.query_alignment_end - record.query_alignment_start)/get_read_length(record.cigarstring)
                if record.mapq < 60 or aligned_frac < 0.7:
                    continue
                read_pos = get_insert_read_pos(record.cigarstring, record.reference_start, pos)
                reads[2][contig+":"+str(pos)].append([record.query_name, read_pos])
            if read_count < 5:
                continue
            insert_contig_pos[contig][pos:pos+1] = 1
            count += 1
        except ValueError:
            pass

# For each read in reads at each pos, randomly select some fraction between 1 and 50% to add insert to.
# If we've seen read already, skip
print("Sample\tRead\tOriginalPosition\tCutFromStart\tCutFromEnd\tInsertLength\tAnnotation\tSourceRead\tSourcePosition\tContigPos\tSeq")
seen_reads = {}
with open(args.output_fastq_1, 'w') as fout_1, open(args.output_fastq_2, 'w') as fout_2:
    with pysam.FastaFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_3_"+args.fastq_suffix) as fin_1:
        for pos in reads[1]:
            # Select an insertion pos
            insert_pos = random.choice(list(inserts.keys()))
            n = random.randint(1, int(0.5*len(reads[1][pos]))+1)
            for i in range(n):
                if reads[1][pos][i][0] in seen_reads:
                    continue
                seq = fin_1.fetch(reads[1][pos][i][0])
                insert = random.choice(inserts[insert_pos])
                insert_seq = insert[5]
                insert_type = insert[3]
                # Truncate the new seq
                new_seq = seq[0:reads[1][pos][i][1]]+insert_seq+seq[reads[1][pos][i][1]:]
                cut = random.randint(1,len(insert_seq))
                new_seq = new_seq[cut:]
                new_seq = new_seq[:len(new_seq)-(len(insert_seq)-cut)]
                qual = '='*len(new_seq)
                fout_1.write("@"+reads[1][pos][i][0]+"\n")
                fout_1.write(new_seq+"\n+\n")
                fout_1.write(qual+"\n")
                seen_reads[reads[1][pos][i][0]] = 1
                print("3\t"+reads[1][pos][i][0]+"\t"+str(reads[1][pos][i][1])+"\t"+str(cut)+"\t"+str(len(insert_seq)-cut)+"\t"+str(len(insert_seq))+"\t"+insert_type+"\t"+insert[0]+"\t"+insert[2]+"\t"+pos+"\t"+insert_seq)
    with pysam.FastaFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_5_"+args.fastq_suffix) as fin_2:
        for pos in reads[2]:
            # Select an insertion pos
            insert_pos = random.choice(list(inserts.keys()))
            n = random.randint(1, int(0.5*len(reads[1][pos]))+1)
            for i in range(n):
                if reads[2][pos][i][0] in seen_reads:
                    continue
                seq = fin_2.fetch(reads[2][pos][i][0])
                insert = random.choice(inserts[insert_pos])
                insert_seq = insert[5]
                insert_type = insert[3]
                new_seq = seq[0:reads[2][pos][i][1]]+insert_seq+seq[reads[2][pos][i][1]:]
                cut = random.randint(1,len(insert_seq))
                new_seq = new_seq[cut:]
                new_seq = new_seq[:len(new_seq)-(len(insert_seq)-cut)]
                qual = '='*len(new_seq)
                fout_2.write("@"+reads[2][pos][i][0]+"\n")
                fout_2.write(new_seq+"\n+\n")
                fout_2.write(qual+"\n")
                seen_reads[reads[2][pos][i][0]] = 1
                print("5\t"+reads[2][pos][i][0]+"\t"+str(reads[2][pos][i][1])+"\t"+str(cut)+"\t"+str(len(insert_seq)-cut)+"\t"+str(len(insert_seq))+"\t"+insert_type+"\t"+insert[0]+"\t"+insert[2]+"\t"+pos+"\t"+insert_seq)
    # Once we have the reads for each sample open up the files and begin writing
    with pysam.FastxFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_3_"+args.fastq_suffix) as fin_1:
        for entry in fin_1:
            if entry.name not in seen_reads:
                fout_1.write("@"+entry.name+"\n")
                fout_1.write(entry.sequence+"\n+\n")
                fout_1.write(entry.quality+"\n")
                
    with pysam.FastxFile(args.sample+"/"+args.fastq_folder+"/"+args.sample+"_5_"+args.fastq_suffix) as fin_2:
        for entry in fin_2:
            if entry.name not in seen_reads:
                fout_2.write("@"+entry.name+"\n")
                fout_2.write(entry.sequence+"\n+\n")
                fout_2.write(entry.quality+"\n")

