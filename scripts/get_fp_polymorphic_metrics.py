import argparse
import sys
import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import pysam
from pysam import VariantFile

parser = argparse.ArgumentParser( description='Get metrics for runs')
parser.add_argument('--input-annotation', required=True)
parser.add_argument('--input-prefix', required=True)
parser.add_argument('--id', required=True)
parser.add_argument('--folder', required=True)
parser.add_argument('--suffix', default=".all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv")
parser.add_argument('--filtered-suffix', default=".filter_tsvs.tsv")
parser.add_argument('--sv-folder', default="structural_variants")
parser.add_argument('--sniffles-suffix', default=".sniffles.vcf")
parser.add_argument('--cutesv-suffix', default=".cuteSV.vcf")
parser.add_argument('--sniffles-filtered-suffix', default=".sniffles.annotated.vcf")
parser.add_argument('--cutesv-filtered-suffix', default=".cuteSV.annotated.vcf")

args = parser.parse_args()

### Read Level ###

# Open up the annotation file
truth_set = defaultdict(list)
truth_positions = defaultdict(IntervalTree)
count = 0
with open(args.input_annotation, 'r') as in_an:
    for line in in_an:
        if count == 0:
            count = 1
            continue
        row = line.strip().split('\t')
        truth_set[row[0]] = row
        chrom = row[2].split(':')[0]
        start = int(row[2].split(':')[1].split('-')[0])
        end = int(row[2].split(':')[1].split('-')[1])
        truth_positions[chrom][(start-500):(end+500)] = row[3]+"_"+row[4]

# Check Pipeline first
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(list)
tn = 0

for sample in ["Sample3", "Sample5"]:
    with open(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.folder+"/"+sample+args.suffix, 'r') as in_data:
        for line in in_data:
            row = line.strip().split('\t')
            if row[3] in truth_set:
                # have a read that we care about
                if row[8] == "PASS":
                    # Check more here
                    tp[row[3]] = 1
                else:
                    fn[row[3]].append(line.strip())
            else:
                # Read we don't care about
                if row[8] == "PASS":
                    fp[row[3]] += 1
                else:
                    tn += 1

# now check to see if we've missed anything.
# Look for reads in truth set not found in 
missed_fn = defaultdict(int)
for read in truth_set:
    if read in tp:
        continue
    if read not in fn:
        missed_fn[read] = 1
        #print("Missed\t"+"\t".join(truth_set[read]))

# Remove any FN calls if we already have seen a TP for that read
# Possible the read has more than one insertion and we have found the one we care abot
# Other calls could be polymorphic or not RT inserts being called as FNs
fn_count = 0
for read in fn:
    if read not in tp:
        fn_count += 1
        #print("\t".join(truth_set[read]))
        #for item in fn[read]:
            #print("found\t"+item)

print("Pipeline")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")



# Check Pipeline with no filters
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(list)
tn = 0

for sample in ["Sample3", "Sample5"]:
    with open(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.folder+"/"+sample+args.filtered_suffix, 'r') as in_data:
        for line in in_data:
            row = line.strip().split('\t')
            if row[3] in truth_set:
                # have a read that we care about
                if row[8] == "PASS":
                    # Check more here
                    tp[row[3]] = 1
                else:
                    fn[row[3]].append(line.strip())
            else:
                # Read we don't care about
                if row[8] == "PASS":
                    fp[row[3]] += 1
                else:
                    tn += 1

# now check to see if we've missed anything.
# Look for reads in truth set not found in 
missed_fn = defaultdict(int)
for read in truth_set:
    if read in tp:
        continue
    if read not in fn:
        missed_fn[read] = 1
        #print("Missed\t"+"\t".join(truth_set[read]))

# Remove any FN calls if we already have seen a TP for that read
# Possible the read has more than one insertion and we have found the one we care abot
# Other calls could be polymorphic or not RT inserts being called as FNs
fn_count = 0
for read in fn:
    if read not in tp:
        fn_count += 1
        #print("\t".join(truth_set[read]))
        #for item in fn[read]:
            #print("found\t"+item)

print("Pipeline")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")



# Now check Sniffles
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(int)
tn = 0

# Load up Sample 1 first and then use it to filter within 500bp with any insertion >=50 bp in Sample3/5
# Need to do this to mimic use of Sample 1 for control filtering
control_inserts_per_chrom = defaultdict(IntervalTree)
sniffles_control = pysam.VariantFile(args.input_prefix+"/"+args.id+"/Sample1/"+args.sv_folder+"/Sample1"+args.sniffles_suffix)
for rec in sniffles_control.fetch():
    if rec.info["SVTYPE"] == "INS" and int(rec.info["SVLEN"]) >= 50:
        # Add this to the dict
        control_inserts_per_chrom[rec.chrom][int(rec.start-500):int(rec.stop+500)] = 1

for sample in ["Sample3", "Sample5"]:
    sniffles_sample = pysam.VariantFile(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.sv_folder+"/"+sample+args.sniffles_filtered_suffix)
    for rec in sniffles_sample.fetch():
        if int(rec.info["SVLEN"]) >= 100:
            insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
            if len(insert_list) == 0:
                # Insert does not appear in Sample 1
                for name in rec.info["RNAMES"]:
                    if name in truth_set:
                        # Have a true positive here
                        tp[name] = 1
                    else:
                        fp[name] = 1
            else:
                # Insertion expected but found to be polymorphic in Sample 1
                for name in rec.info["RNAMES"]:
                    if name in truth_set:
                        fn[name] = 1
                    else:
                        tn += 1

missed_fn = defaultdict(int)
for read in truth_set:
    if read in tp:
        continue
    if read not in fn:
        missed_fn[read] = 1

fn_count = 0
for read in fn:
    if read not in tp:
        fn_count += 1

print("Sniffles")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")


# Now check CuteSV
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(int)
tn = 0

# Load up Sample 1 first and then use it to filter within 500bp with any insertion >=50 bp in Sample3/5
# Need to do this to mimic use of Sample 1 for control filtering
control_inserts_per_chrom = defaultdict(IntervalTree)
cuteSV_control = pysam.VariantFile(args.input_prefix+"/"+args.id+"/Sample1/"+args.sv_folder+"/Sample1"+args.cutesv_suffix)
for rec in cuteSV_control.fetch():
    if rec.info["SVTYPE"] == "INS" and int(rec.info["SVLEN"]) >= 50:
        # Add this to the dict
        control_inserts_per_chrom[rec.chrom][int(rec.start-500):int(rec.stop+500)] = 1

for sample in ["Sample3", "Sample5"]:
    cuteSV_sample = pysam.VariantFile(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.sv_folder+"/"+sample+args.cutesv_filtered_suffix)
    for rec in cuteSV_sample.fetch():
        if int(rec.info["SVLEN"]) >= 100:
            insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
            if len(insert_list) == 0:
                # Insert does not appear in Sample 1
                for name in rec.info["RNAMES"]:
                    if name in truth_set:
                        # Have a true positive here
                        tp[name] = 1
                    else:
                        fp[name] = 1
            else:
                # Insertion expected but found to be polymorphic in Sample 1
                for name in rec.info["RNAMES"]:
                    if name in truth_set:
                        fn[name] = 1
                    else:
                        tn += 1

missed_fn = defaultdict(int)
for read in truth_set:
    if read in tp:
        continue
    if read not in fn:
        missed_fn[read] = 1

fn_count = 0
for read in fn:
    if read not in tp:
        fn_count += 1

print("CuteSV")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")

### Position Level ###

# Check Pipeline first
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(list)
tn = 0

for sample in ["Sample3", "Sample5"]:
    with open(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.folder+"/"+sample+args.suffix, 'r') as in_data:
        count = 0
        for line in in_data:
            if count == 0:
                count = 1
                continue
            row = line.strip().split('\t')
            nearby = truth_positions[row[0]][int(row[1]):int(row[2])]
            if len(nearby) > 0:
                # have a position that we care about
                if row[8] == "PASS":
                    # Check more here
                    for item in nearby:
                        tp[row[0]+":"+str(item.begin)+"-"+str(item.end)] += 1
                else:
                    for item in nearby:
                        fn[row[0]+":"+str(item.begin)+"-"+str(item.end)].append(line.strip())
            else:
                # Read we don't care about
                if row[8] == "PASS":
                    fp[row[0]+":"+row[1]+"-i"+row[2]] += 1
                else:
                    tn += 1

# now check to see if we've missed anything.
# Look for reads in truth set not found in 
missed_fn = defaultdict(int)
for chrom in truth_positions:
    for item in truth_positions[chrom]:
        if chrom+":"+str(item.begin)+"-"+str(item.end) in tp:
            continue
        if chrom+":"+str(item.begin)+"-"+str(item.end) not in fn:
            missed_fn[chrom+":"+str(item.begin)+"-"+str(item.end)] = 1

# Remove any FN calls if we already have seen a TP for that read
# Possible the read has more than one insertion and we have found the one we care abot
# Other calls could be polymorphic or not RT inserts being called as FNs
fn_count = 0
for pos in fn:
    if pos not in tp:
        fn_count += 1
        #print("\t".join(truth_set[read]))
        #for item in fn[read]:
            #print("found\t"+item)

print("Position Level Pipeline")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")


# Check Pipeline with no filters
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(list)
tn = 0

for sample in ["Sample3", "Sample5"]:
    with open(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.folder+"/"+sample+args.filtered_suffix, 'r') as in_data:
        count = 0
        for line in in_data:
            if count == 0:
                count = 1
                continue
            row = line.strip().split('\t')
            nearby = truth_positions[row[0]][int(row[1]):int(row[2])]
            if len(nearby) > 0:
                # have a position that we care about
                if row[8] == "PASS":
                    # Check more here
                    for item in nearby:
                        tp[row[0]+":"+str(item.begin)+"-"+str(item.end)] += 1
                else:
                    for item in nearby:
                        fn[row[0]+":"+str(item.begin)+"-"+str(item.end)].append(line.strip())
            else:
                # Read we don't care about
                if row[8] == "PASS":
                    fp[row[0]+":"+row[1]+"-i"+row[2]] += 1
                else:
                    tn += 1

# now check to see if we've missed anything.
# Look for reads in truth set not found in 
missed_fn = defaultdict(int)
for chrom in truth_positions:
    for item in truth_positions[chrom]:
        if chrom+":"+str(item.begin)+"-"+str(item.end) in tp:
            continue
        if chrom+":"+str(item.begin)+"-"+str(item.end) not in fn:
            missed_fn[chrom+":"+str(item.begin)+"-"+str(item.end)] = 1

# Remove any FN calls if we already have seen a TP for that read
# Possible the read has more than one insertion and we have found the one we care abot
# Other calls could be polymorphic or not RT inserts being called as FNs
fn_count = 0
for pos in fn:
    if pos not in tp:
        fn_count += 1
        #print("\t".join(truth_set[read]))
        #for item in fn[read]:
            #print("found\t"+item)

print("Position Level Pipeline")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")


# Now check Sniffles
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(int)
tn = 0

# Load up Sample 1 first and then use it to filter within 500bp with any insertion >=50 bp in Sample3/5
# Need to do this to mimic use of Sample 1 for control filtering
control_inserts_per_chrom = defaultdict(IntervalTree)
sniffles_control = pysam.VariantFile(args.input_prefix+"/"+args.id+"/Sample1/"+args.sv_folder+"/Sample1"+args.sniffles_suffix)
for rec in sniffles_control.fetch():
    if rec.info["SVTYPE"] == "INS" and int(rec.info["SVLEN"]) >= 50:
        # Add this to the dict
        control_inserts_per_chrom[rec.chrom][int(rec.start-500):int(rec.stop+500)] = 1

for sample in ["Sample3", "Sample5"]:
    sniffles_sample = pysam.VariantFile(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.sv_folder+"/"+sample+args.sniffles_filtered_suffix)
    for rec in sniffles_sample.fetch():
        if int(rec.info["SVLEN"]) >= 100:
            insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
            if len(insert_list) == 0:
                # Insert does not appear in Sample 1
                nearby = truth_positions[rec.chrom][rec.start:rec.stop]
                if len(nearby) > 0:
                    for item in nearby:
                        tp[row[0]+":"+str(item.begin)+"-"+str(item.end)] += 1
                else:
                    fp[rec.chrom+":"+str(rec.start)+"-"+str(rec.stop)] = 1
            else:
                # Insertion expected but found to be polymorphic in Sample 1
                nearby = truth_positions[rec.chrom][rec.start:rec.stop]
                if len(nearby) > 0:
                    for item in nearby:
                        fn[row[0]+":"+str(item.begin)+"-"+str(item.end)] += 1
                else:
                    tn += 1

missed_fn = defaultdict(int)
for chrom in truth_positions:
    for item in truth_positions[chrom]:
        if chrom+":"+str(item.begin)+"-"+str(item.end) in tp:
            continue
        if chrom+":"+str(item.begin)+"-"+str(item.end) not in fn:
            missed_fn[chrom+":"+str(item.begin)+"-"+str(item.end)] = 1

fn_count = 0
for pos in fn:
    if pos not in tp:
        fn_count += 1


print("Position Level Sniffles")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")


# Now check CuteSV
tp = defaultdict(int)
fp = defaultdict(int)
fn = defaultdict(int)
tn = 0

# Load up Sample 1 first and then use it to filter within 500bp with any insertion >=50 bp in Sample3/5
# Need to do this to mimic use of Sample 1 for control filtering
control_inserts_per_chrom = defaultdict(IntervalTree)
cuteSV_control = pysam.VariantFile(args.input_prefix+"/"+args.id+"/Sample1/"+args.sv_folder+"/Sample1"+args.cutesv_suffix)
for rec in cuteSV_control.fetch():
    if rec.info["SVTYPE"] == "INS" and int(rec.info["SVLEN"]) >= 50:
        # Add this to the dict
        control_inserts_per_chrom[rec.chrom][int(rec.start-500):int(rec.stop+500)] = 1

for sample in ["Sample3", "Sample5"]:
    cuteSV_sample = pysam.VariantFile(args.input_prefix+"/"+args.id+"/"+sample+"/"+args.sv_folder+"/"+sample+args.cutesv_filtered_suffix)
    for rec in cuteSV_sample.fetch():
        if int(rec.info["SVLEN"]) >= 100:
            insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
            if len(insert_list) == 0:
                # Insert does not appear in Sample 1
                nearby = truth_positions[rec.chrom][rec.start:rec.stop]
                if len(nearby) > 0:
                    for item in nearby:
                        tp[row[0]+":"+str(item.begin)+"-"+str(item.end)] += 1
                else:
                    fp[rec.chrom+":"+str(rec.start)+"-"+str(rec.stop)] = 1
            else:
                # Insertion expected but found to be polymorphic in Sample 1
                nearby = truth_positions[rec.chrom][rec.start:rec.stop]
                if len(nearby) > 0:
                    for item in nearby:
                        fn[row[0]+":"+str(item.begin)+"-"+str(item.end)] += 1
                else:
                    tn += 1

missed_fn = defaultdict(int)
for chrom in truth_positions:
    for item in truth_positions[chrom]:
        if chrom+":"+str(item.begin)+"-"+str(item.end) in tp:
            continue
        if chrom+":"+str(item.begin)+"-"+str(item.end) not in fn:
            missed_fn[chrom+":"+str(item.begin)+"-"+str(item.end)] = 1

fn_count = 0
for pos in fn:
    if pos not in tp:
        fn_count += 1

print("Position Level CuteSV")
print("TP: "+str(len(tp)))
print("FP: "+str(len(fp)))
print("FN: "+str(fn_count))
print("Missed FN: "+str(len(missed_fn)))
print("Total FN: "+str(fn_count+len(missed_fn)))
print("TN: "+str(tn)+"\n")