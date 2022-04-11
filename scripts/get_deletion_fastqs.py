import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
import gzip
import random
from intervaltree import Interval, IntervalTree

random.seed(10)

def get_read_length(cigarstring):
    # Gets the read length based on the Cigar String
    count = 0
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
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
    return count


def get_insert_position_for_read(record, ref_start, ref_end):
    ref_count = record.reference_start
    prev_ref_count = 0
    failed_start = False
    read_count = 0
    read_start = -1
    read_end = -1
    for cg in re.findall('[0-9]*[A-Z=]', record.cigarstring):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith('D'):
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        if ref_count > ref_start and read_start < 0:
            read_start = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_start = read_start - (ref_count - ref_start)
            if prev_ref_count < 500:
                failed_start = True
        if ref_count > ref_end and read_end < 0:
            read_end = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_end = read_end - (ref_count - ref_end)
    if ref_count - ref_end < 500 or failed_start:
        read_end = -1
    return (read_start, read_end)



def get_insert_position_for_control(record, ref_start, ref_end):
    ref_count = record.reference_start
    prev_ref_count = 0
    failed_start = False
    read_count = 0
    read_start = -1
    read_end = -1
    positions = {}
    count = 0
    for cg in re.findall('[0-9]*[A-Z=]', record.cigarstring):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            start = -1
            end = -1
            if(int(cg[:cg.find("I")])) >= 100 and ref_count > ref_start and ref_count < ref_end:
                start = read_count
            read_count += int(cg[:cg.find("I")])
            if(int(cg[:cg.find("I")])) >= 100 and ref_count > ref_start and ref_count < ref_end:
                end = read_count
                positions[count] = [start,end]
                count += 1
        elif cg.endswith('D'):
            prev_ref_count = ref_count
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
        if ref_count > ref_start and read_start < 0:
            read_start = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_start = read_start - (ref_count - ref_start)
            if prev_ref_count < 500:
                failed_start = True
        if ref_count > ref_end and read_end < 0:
            read_end = read_count
            if cg.endswith('M') or cg.endswith('X') or cg.endswith('='):
                read_end = read_end - (ref_count - ref_end)
    if ref_count - ref_end < 500 or failed_start:
        read_end = -1
    return positions

def mapped_end_to_end(record, read_length):
    ref_count = record.reference_start
    read_count = 0
    read_start = -1
    read_end = -1
    large_indel = False
    for cg in re.findall('[0-9]*[A-Z=]', record.cigarstring):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
        elif cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
            if int(cg[:cg.find("I")]) >= 100:
                large_indel = True
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
            if int(cg[:cg.find("D")]) >= 100:
                large_indel = True
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
            if read_start == -1:
                read_start = read_count
            elif read_end == -1:
                read_end = read_count - int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
            if read_start == -1:
                read_start = read_count
            elif read_end == -1:
                read_end = read_count - int(cg[:cg.find("H")])
        if read_start == -1:
            read_start = read_count
    if read_end == -1:
        read_end = read_count
    if (read_end - read_start)/read_length < 0.5:
        # Need at least half the read to map
        return False
    if large_indel:
        return False
    return True


def get_inserts_for_control(ref_bam, mat_and_pat_inserts, chr_list):
    ret = {}
    bam = pysam.AlignmentFile(ref_bam, "rb")
    for chrom in mat_and_pat_inserts:
        if chrom in chr_list:
            for item in mat_and_pat_inserts[chrom]:
                insert_dict = item.data
                read_positions = {}
                to_print = False
                if "mat" in insert_dict:
                    mat_pos = insert_dict["mat"][0]+":"+str(int(insert_dict["mat"][1])-500)+"-"+str(int(insert_dict["mat"][2])+500)
                    for record in bam.fetch(region=mat_pos):
                        # Check read mapped end to end with no inserts
                        positions = get_insert_position_for_control(record, int(insert_dict["mat"][1])-500, int(insert_dict["mat"][2])+500)
                        if len(positions) == 0:
                            continue
                        read_start = -1
                        read_end = -1
                        for p in positions:
                            if read_start == -1:
                                read_start = positions[p][0]
                            if read_end == -1:
                                read_end = positions[p][1]
                            if positions[p][0] < read_start:
                                read_start = positions[p][0]
                            if positions[p][1] > read_end:
                                read_end = positions[p][1]
                        read_length = get_read_length(record.cigarstring)
                        if read_start < 0 or read_end < 0:
                            continue
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        read_positions[record.query_name] = [read_start, read_end, record.cigarstring, record.is_reverse]
                if "pat" in insert_dict:
                    pat_pos = insert_dict["pat"][0]+":"+str(int(insert_dict["pat"][1])-500)+"-"+str(int(insert_dict["pat"][2])+500)
                    for record in bam.fetch(region=pat_pos):
                        positions = get_insert_position_for_control(record, int(insert_dict["pat"][1])-500, int(insert_dict["pat"][2])+500)
                        if len(positions) == 0:
                            continue
                        read_start = -1
                        read_end = -1
                        for p in positions:
                            if read_start == -1:
                                read_start = positions[p][0]
                            if read_end == -1:
                                read_end = positions[p][1]
                            if positions[p][0] < read_start:
                                read_start = positions[p][0]
                            if positions[p][1] > read_end:
                                read_end = positions[p][1]
                        read_length = get_read_length(record.cigarstring)
                        if read_start < 0 or read_end < 0:
                            continue
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        read_positions[record.query_name] = [read_start, read_end, record.cigarstring, record.is_reverse]
                ret[chrom+":"+str(item.begin)+"-"+str(item.end)] = read_positions
    return ret



def get_inserts_for_sample(mat_bam, pat_bam, mat_and_pat_inserts, chr_list, file):
    ret = {}
    depths = {}
    mat = pysam.AlignmentFile(mat_bam, "rb")
    pat = pysam.AlignmentFile(pat_bam, "rb")
    for chrom in mat_and_pat_inserts:
        if chrom in chr_list:
            for item in mat_and_pat_inserts[chrom]:
                insert_dict = item.data
                read_positions = {}
                to_print = False
                #if chrom == "chr7" and item.begin > 85262406 and item.end < 85282407:
                #    to_print = True
                depth = 0
                if "mat" in insert_dict and "pat" in insert_dict:
                    # Get the reads that aligned to the mat first, and then for pat
                    mat_pos = insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4])-100)+"-"+str(int(insert_dict["mat"][5])+100)
                    for record in mat.fetch(region=mat_pos):
                        read_length = get_read_length(record.cigarstring)
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["mat"][4]), int(insert_dict["mat"][5]))
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        if read_start < 0 or read_end < 0:
                            continue
                        depth += 1
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        #if record.query_name == "90f746b1-3df1-4d99-989c-4507e3a3ca28":
                        #    print(record.query_name+"\t"+str(read_start)+"\t"+str(read_end))
                        read_positions[record.query_name] = [read_start, read_end, insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4]))+"-"+str(int(insert_dict["mat"][5]))]
                    pat_pos = insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4])-100)+"-"+str(int(insert_dict["pat"][5])+100)
                    for record in pat.fetch(region=pat_pos):
                        read_length = get_read_length(record.cigarstring)
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["pat"][4]), int(insert_dict["pat"][5]))
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        if read_start < 0 or read_end < 0:
                            continue
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        if record.query_name not in read_positions:
                            depth += 1
                            #if record.query_name == "90f746b1-3df1-4d99-989c-4507e3a3ca28":
                            #   print(record.query_name+"\t"+str(read_start)+"\t"+str(read_end))
                            read_positions[record.query_name] = [read_start, read_end, insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4]))+"-"+str(int(insert_dict["pat"][5]))]
                elif "mat" in insert_dict:
                    mat_pos = insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4])-100)+"-"+str(int(insert_dict["mat"][5])+100)
                    for record in mat.fetch(region=mat_pos):
                        # Check read mapped end to end with no inserts
                        read_length = get_read_length(record.cigarstring)
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        #if not mapped_end_to_end(record, read_length):
                        #    continue
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["mat"][4]), int(insert_dict["mat"][5]))
                        if read_start < 0 or read_end < 0:
                            continue
                        depth += 1
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        #if record.query_name == "90f746b1-3df1-4d99-989c-4507e3a3ca28":
                        #    print(record.query_name+"\t"+str(read_start)+"\t"+str(read_end))
                        read_positions[record.query_name] = [read_start, read_end, insert_dict["mat"][3]+":"+str(int(insert_dict["mat"][4]))+"-"+str(int(insert_dict["mat"][5]))]
                elif "pat" in insert_dict:
                    pat_pos = insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4])-100)+"-"+str(int(insert_dict["pat"][5])+100)
                    for record in pat.fetch(region=pat_pos):
                        read_length = get_read_length(record.cigarstring)
                        if file is not None:
                            with open(file, "a") as myfile:
                                myfile.write(record.query_name+"\n")
                        #if not mapped_end_to_end(record, read_length):
                        #    continue
                        read_start, read_end = get_insert_position_for_read(record, int(insert_dict["pat"][4]), int(insert_dict["pat"][5]))
                        if read_start < 0 or read_end < 0:
                            continue
                        depth += 1
                        if record.is_reverse:
                            read_start = read_length-read_start
                            read_end = read_length - read_end
                        if read_start > read_end:
                            tmp = read_start
                            read_start = read_end
                            read_end = tmp
                        read_positions[record.query_name] = [read_start, read_end, insert_dict["pat"][3]+":"+str(int(insert_dict["pat"][4]))+"-"+str(int(insert_dict["pat"][5]))]
                ret[chrom+":"+str(item.begin)+"-"+str(item.end)] = read_positions
                depths[chrom+":"+str(item.begin)+"-"+str(item.end)] = depth
    return ret, depths


def select_reads_to_delete(insert_positions, delete_all, deletion_fraction, depth_per_pos, novel_positions, sample):
    read_positions_to_delete = defaultdict(list)
    read_positions_to_retain = defaultdict(list)
    count_and_depth_per_pos = defaultdict(str)
    for pos in insert_positions:
        to_print = False
        #if pos.split(':')[0] == "chr18" and int(pos.split(':')[1].split('-')[0]) > 40700659 and int(pos.split(':')[1].split('-')[1]) < 40768732:
        #    to_print = True
         #For each position either delete all or some fraction
        count = 0
        novel = False
        s = ""
        if pos in novel_positions:
            novel = True
            s = novel_positions[pos]
        for read in insert_positions[pos]:
            if delete_all:
                read_positions_to_delete[read].append(insert_positions[pos][read])
                if to_print:
                    print("Deleted: "+"\t"+read+"\t"+str(insert_positions[pos][read][0])+"\t"+str(insert_positions[pos][read][1]))
                    print(insert_positions[pos][read][2])
                    print(insert_positions[pos][read][3])
            elif novel:
                # Only keep one read for the position
                read_start = insert_positions[pos][read][0]
                read_end = insert_positions[pos][read][1]
                if count == 0 and s == sample:
                    read_positions_to_retain[read].append([read_start, read_end, pos])
                    if to_print:
                        print("Retained: "+"\t"+read+"\t"+str(read_start)+"\t"+str(read_end))
                    count += 1
                else:
                    if to_print:
                        print("Deleted: "+"\t"+read+"\t"+str(read_start)+"\t"+str(read_end))
                    read_positions_to_delete[read].append([read_start, read_end])
            else:
                # Have to randomly select which reads to delete
                read_start = insert_positions[pos][read][0]
                read_end = insert_positions[pos][read][1]
                if count < 1:
                    # Have at least two read per position
                    read_positions_to_retain[read].append([read_start, read_end, pos])
                    if to_print:
                        print("Retained: "+"\t"+read+"\t"+str(read_start)+"\t"+str(read_end))
                    count += 1
                elif random.uniform(0, 1) < deletion_fraction:
                    if to_print:
                        print("Deleted: "+"\t"+read+"\t"+str(read_start)+"\t"+str(read_end))
                    read_positions_to_delete[read].append([read_start, read_end])
                else:
                    # These are the reads and positions we need to retain
                    read_positions_to_retain[read].append([read_start, read_end, pos])
                    if to_print:
                        print("Retained: "+"\t"+read+"\t"+str(read_start)+"\t"+str(read_end))
        if not delete_all:
            count_and_depth_per_pos[pos] = str(count)+"\t"+str(depth_per_pos[pos])
    return (read_positions_to_delete, read_positions_to_retain, count_and_depth_per_pos)

def get_fastq(sample, args):
    if sample == "Sample1":
        return args.input_fastq_1
    if sample == "Sample3":
        return args.input_fastq_3
    if sample == "Sample5":
        return args.input_fastq_5
    return None

def get_output_fastq(sample, args):
    if sample in ["Sample1", "Sample3", "Sample5"]:
        return args.output_fasta_folder+"/"+sample+".fastq"
    return None

parser = argparse.ArgumentParser( description='Generate New fastqs with deletions in reads')
parser.add_argument('--output-fasta-folder', default="fastq")
parser.add_argument('--input-fastq-1', required=True)
parser.add_argument('--input-fastq-3', required=True)
parser.add_argument('--input-fastq-5', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--centromeres', required=True)
parser.add_argument('--chr-list', default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY")
parser.add_argument('--mat-ref-inserts', required=True)
parser.add_argument('--pat-ref-inserts', required=True)
parser.add_argument('--mat-bam-1', required=True)
parser.add_argument('--mat-bam-3', required=True)
parser.add_argument('--mat-bam-5', required=True)
parser.add_argument('--pat-bam-1', required=True)
parser.add_argument('--pat-bam-3', required=True)
parser.add_argument('--pat-bam-5', required=True)
parser.add_argument('--reads-to-ignore', required=True)
args = parser.parse_args()

#print("centromeres")
# Read in Centromeres list
centromeres = defaultdict(IntervalTree)
with open(args.centromeres) as in_cf:
    count = 0
    for line in in_cf:
        if count == 0:
            count = 1
            continue
        line_args = line.strip().split('\t')
        chrom = line_args[1]
        start = int(line_args[2])
        end = int(line_args[3])
        key = chrom + ":" + str(start) + "-" + str(end)
        centromeres[chrom][start:end] = key

#print("mat")
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
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        nearby = centromeres[chrom][start-500:end+500]
        if row[8] != "PASS" or len(nearby) > 0:
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
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        nearby = centromeres[chrom][start-500:end+500]
        if row[8] != "PASS" or len(nearby) > 0:
            continue
        pat_ref_inserts[chrom][start:end] = row
        nearby = mat_ref_inserts[chrom][start-500:end+500]
        if len(nearby) > 0:
            # Add both the mat and pat inserts to the combined at this position, should only have one item given the polished assembly
            for item in nearby:
                mat_row = item.data
                combined_mat_pat_ref_inserts[chrom][start:end] = {"pat": row, "mat": mat_row}
        else:
            # Add just the pat row
            combined_mat_pat_ref_inserts[chrom][start:end] = {"pat": row}

#print("mat2")
count = 0
with open(args.mat_ref_inserts, 'r') as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        nearby = centromeres[chrom][start-500:end+500]
        if row[8] != "PASS" or len(nearby) > 0:
            continue
        nearby = pat_ref_inserts[chrom][start-500:end+500]
        if len(nearby) > 0:
            # Already added this insert to the combined set
            continue
        else:
            combined_mat_pat_ref_inserts[chrom][start:end] = {"mat": row}


# Read the bams for each sample against the mat inserts and the pat inserts. Decide which reads map to each sample at the position. 

# Now read in each of the read fastq inserts
# If the insert passes and lines up within 500 bp of an insert in the paternal contig set, flag it
# In sample 1 instantly collapse these to generate an insert free set, while in sample 3 and 5 collapse some fraction of the reads at each position
reads_with_inserts_per_sample = {}
depth_per_sample = {}
reads_with_inserts_per_sample["Sample1"], depth_per_sample["Sample1"] = get_inserts_for_sample(args.mat_bam_1, args.pat_bam_1, combined_mat_pat_ref_inserts, args.chr_list.split(","),args.reads_to_ignore)
reads_with_inserts_per_sample["Sample3"], depth_per_sample["Sample3"] = get_inserts_for_sample(args.mat_bam_3, args.pat_bam_3, combined_mat_pat_ref_inserts, args.chr_list.split(","),None)
reads_with_inserts_per_sample["Sample5"], depth_per_sample["Sample5"] = get_inserts_for_sample(args.mat_bam_5, args.pat_bam_5, combined_mat_pat_ref_inserts, args.chr_list.split(","),None)

# Get positions where we will have novel insertions
novel_positions = defaultdict(int)
for sample in reads_with_inserts_per_sample:
    for pos in reads_with_inserts_per_sample[sample]:
        if pos in novel_positions:
            continue
        if random.uniform(0, 1) < 0.25:
            # Position selected to be novel
            s = "Sample3"
            if random.uniform(0, 1) < 0.5:
                s = "Sample5"
            novel_positions[pos] = s

# Print header
print("Sample\tRead\tReadStart\tReadEnd\tRefPosition")

# Once we have insertions on reads, we now have to decide at each position how many reads to delete sequence and how many to retain
# For Sample 1 delete all inserts at every position - clean control sample
# For Samples 3 and 5 decide randomly based on some fraction of which to retian
for sample in ["Sample1", "Sample3", "Sample5"]:
    #print(sample)
    delete_all = False
    deletion_fraction = 0.9
    if sample == "Sample1":
        delete_all = True
    deletion_positions, retained_positions, count_and_depth_per_pos = select_reads_to_delete(reads_with_inserts_per_sample[sample], delete_all, deletion_fraction, depth_per_sample[sample], novel_positions, sample)
    with open(get_output_fastq(sample, args), 'w') as fout:
        with pysam.FastxFile(get_fastq(sample, args)) as fin:
            for entry in fin:
                seq = entry.sequence
                if entry.name in deletion_positions and entry.name in retained_positions:
                    deletion_positions_list = deletion_positions[entry.name]
                    retained_positions_list = retained_positions[entry.name]
                    deletion_positions_list.sort()
                    retained_positions_list.sort()
                    deleted_count = 0
                    rp_count = 0
                    for i in range(len(deletion_positions_list)):
                        if rp_count < len(retained_positions_list):
                            if deletion_positions_list[i][0] > retained_positions_list[rp_count][0]:
                                # Print out the retained position at rp_count
                                print(sample+"\t"+entry.name+"\t"+str(retained_positions_list[rp_count][0]-deleted_count)+"\t"+str(retained_positions_list[rp_count][1]-deleted_count)+"\t"+str(retained_positions_list[rp_count][2])+"\t"+count_and_depth_per_pos[str(retained_positions_list[rp_count][2])])
                                rp_count+= 1
                        # Delete the seqeunce
                        start = deletion_positions_list[i][0] - deleted_count
                        end = deletion_positions_list[i][1] - deleted_count
                        seq = seq[0:start]+seq[end:]
                        deleted_count += (end - start)
                    # Output any remaining retained positions that are after the deleted positions
                    if len(seq) < 1:
                        continue
                    for i in range(len(retained_positions_list)):
                        if i < rp_count:
                            continue
                        print(sample+"\t"+entry.name+"\t"+str(retained_positions_list[i][0]-deleted_count)+"\t"+str(retained_positions_list[i][1]-deleted_count)+"\t"+str(retained_positions_list[i][2])+"\t"+count_and_depth_per_pos[str(retained_positions_list[i][2])])
                elif entry.name in deletion_positions:
                    deletion_positions_list = deletion_positions[entry.name]
                    deletion_positions_list.sort()
                    deleted_count = 0
                    for i in range(len(deletion_positions_list)):
                        # Delete the seqeunce
                        start = deletion_positions_list[i][0] - deleted_count
                        end = deletion_positions_list[i][1] - deleted_count
                        seq = seq[0:start]+seq[end:]
                        deleted_count += (end - start)
                elif entry.name in retained_positions:
                    retained_positions_list = retained_positions[entry.name]
                    retained_positions_list.sort()
                    for i in range(len(retained_positions_list)):
                        print(sample+"\t"+entry.name+"\t"+str(retained_positions_list[i][0])+"\t"+str(retained_positions_list[i][1])+"\t"+str(retained_positions_list[i][2])+"\t"+count_and_depth_per_pos[str(retained_positions_list[i][2])])
                fout.write("@"+entry.name+"\n")
                fout.write(seq+"\n")
                fout.write("+\n")
                fout.write("="*len(seq)+"\n")




