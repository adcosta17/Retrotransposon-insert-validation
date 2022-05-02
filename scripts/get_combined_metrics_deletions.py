import argparse
import sys
import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import pysam
from pysam import VariantFile

def get_soft_clip(cigartuples, reverse=False):
    # Gets the read length based on the Cigar String
    count = 0
    if reverse:
        cigartuples.reverse()
    for cg in cigartuples:
        if cg[0] == 4:
            count += cg[1]
        elif cg[0] == 5:
            count += cg[1]
        else:
            break
    return count


def get_depth(input_prefix, genome_id, bam_folder, pos, insert_size, samples):
    chrom = pos.split(':')[0]
    start = str(int(pos.split(':')[1].split('-')[0]) - 50)
    end = str(int(pos.split(':')[1].split('-')[1]) + 50)
    hp1 = 0
    hp2 = 0
    hp_all = 0
    hp1_soft = 0
    hp2_soft = 0
    hp_all_soft = 0
    for sample in samples:
        sam_reader = pysam.AlignmentFile(input_prefix+"/"+genome_id+"/"+sample+"/"+bam_folder+"/"+sample+".sorted.phased.bam")
        tmp_sam_reader = sam_reader.fetch(region=chrom+":"+start+'-'+end)
        for record in tmp_sam_reader:
            hp_all += 1
            tags = record.get_tags()
            hp = 0
            for t in tags:
                if t[0] == "HP":
                    hp = int(t[1])
            if hp == 0:
                hp1 += 1
                hp2 += 1
            elif hp == 1:
                hp1 += 1
            elif hp == 2:
                hp2 += 1
            if abs(record.reference_start - int(end)) > abs(record.reference_end - int(start)):
                # Look at the end of the CIGAR
                size = get_soft_clip(record.cigartuples, True)
                if size > 0.25*insert_size:
                    hp_all_soft += 1
                    if hp == 0:
                        hp1_soft += 1
                        hp2_soft += 1
                    elif hp == 1:
                        hp1_soft += 1
                    elif hp == 2:
                        hp2_soft += 1
            else:
                # Look at the start
                size = get_soft_clip(record.cigartuples)
                if size > 0.25*insert_size:
                    hp_all_soft += 1
                    if hp == 0:
                        hp1_soft += 1
                        hp2_soft += 1
                    elif hp == 1:
                        hp1_soft += 1
                    elif hp == 2:
                        hp2_soft += 1
    return hp1, hp2, hp_all, hp1_soft, hp2_soft, hp_all_soft

parser = argparse.ArgumentParser( description='Get metrics for runs')
parser.add_argument('--input-prefix', required=True)
parser.add_argument('--id-list', required=True)
parser.add_argument('--folder', required=True)
parser.add_argument('--bam-folder', required=True)
parser.add_argument('--annotation-suffix', default=".annotation.tsv")
parser.add_argument('--annotation-folder', default="fastq_with_deletions")
parser.add_argument('--suffix', default=".all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv")
parser.add_argument('--samples', default="Sample3,Sample5")
parser.add_argument('--sv-folder', default="structural_variants")
parser.add_argument('--sniffles-suffix', default=".sniffles.vcf")
parser.add_argument('--cutesv-suffix', default=".cuteSV.vcf")
parser.add_argument('--sniffles-filtered-suffix', default=".sniffles.annotated.vcf")
parser.add_argument('--cutesv-filtered-suffix', default=".cuteSV.annotated.vcf")
parser.add_argument('--fn', required=False)

args = parser.parse_args()

### Read Level ###
tp_overall = 0
tp_insert_size = defaultdict(int)
tp_read_size = defaultdict(int)
tp_repeat_type = defaultdict(int)
tp_spliced_depth = defaultdict(int)
fp_insert_size = defaultdict(int)
fp_overall = 0
fp_read_size = defaultdict(int)
fp_repeat_type = defaultdict(int)
fp_spliced_depth = defaultdict(int)
fn_insert_size = defaultdict(int)
fn_read_size = defaultdict(int)
fn_repeat_type = defaultdict(int)
fn_spliced_depth = defaultdict(int)
fn_overall = 0
missed_fn_overall = 0
tn = {}

tp_position_overall = 0
tp_position_insert_size = defaultdict(int)
tp_position_repeat_type = defaultdict(int)
tp_position_spliced_depth = defaultdict(int)
fp_position_insert_size = defaultdict(int)
fp_position_overall = 0
fp_position_repeat_type = defaultdict(int)
fp_position_spliced_depth = defaultdict(int)
fn_position_insert_size = defaultdict(int)
fn_position_repeat_type = defaultdict(int)
fn_position_spliced_depth = defaultdict(int)
fn_position_overall = 0
tn_position = 0

fn_rows = defaultdict(list)

for genome_id in args.id_list.split(','):
    tmp = genome_id
    #print(genome_id)
    reads_to_check = defaultdict(int)
    tp = defaultdict(str)
    fp = defaultdict(str)
    fn = defaultdict(str)
    tp_position = defaultdict(str)
    fp_position = defaultdict(str)
    fn_position = defaultdict(str)

    # Get the set of mat and pat insertions. If there is more than 1 insert within 1000bp ignore those positions in the truth data
    mat_inserts = defaultdict(IntervalTree)
    pat_inserts = defaultdict(IntervalTree)
    all_inserts = defaultdict(IntervalTree)
    count = 0
    with open(tmp+"/assembly_subset/"+tmp+".mat.insertions.repbase_annotated.tsv", 'r') as in_mat:
        for line in in_mat:
            if count == 0:
                count =1 
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1]) - 1000
            end = int(row[2]) + 1000
            #mat_all_inserts[chrom][start:end] = 1
            if "PASS" not in line:
                continue
            mat_inserts[chrom][start:end] = 1
    count = 0
    with open(tmp+"/assembly_subset/"+tmp+".pat.insertions.repbase_annotated.tsv", 'r') as in_pat:
        for line in in_pat:
            if count == 0:
                count =1 
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1]) - 1000
            end = int(row[2]) + 1000
            #pat_all_inserts[chrom][start:end] = 1
            if "PASS" not in line:
                continue
            pat_inserts[chrom][start:end] = 1
    count = 0
    with open(tmp+"/assembly_analysis/"+tmp+".insertions.repbase_annotated.tsv", 'r') as in_all:
        for line in in_all:
            if count == 0:
                count =1 
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1]) - 1000
            end = int(row[2]) + 1000
            all_inserts[chrom][start:end] = 1

    # Open up the annotation file
    truth_set = {}
    for sample in args.samples.split(','):
        truth_set[sample] = defaultdict(list)
    truth_set["Novel"] = defaultdict(list)
    truth_set["HotSpot"] = defaultdict(list)
    truth_positions = defaultdict(IntervalTree)
    reads_to_pos = defaultdict(list)
    read_depth = defaultdict(int)
    count = 0
    cut_start = 0
    cut_end = 0
    depth_per_pos = defaultdict(list)
    count_per_pos = defaultdict(int)
    with open(genome_id+"/"+args.annotation_folder+"/"+genome_id+args.annotation_suffix, 'r') as in_an:
        for line in in_an:
            if count == 0:
                count = 1
                continue
            row = line.strip().split('\t')
            sample = row[0]
            chrom = row[4].split(':')[0]
            start = int(row[4].split(':')[1].split('-')[0])
            end = int(row[4].split(':')[1].split('-')[1])
            if len(mat_inserts[chrom][start:end]) > 1 or len(pat_inserts[chrom][start:end]) > 1:
                # Have multiple inserts too close together
                continue
            truth_set[sample][row[1]] = row
            count_per_pos[row[4]] += int(row[5])
            if row[4] not in depth_per_pos:
                hp1, hp2, overall, hp1_soft, hp2_soft, hp_all_soft = get_depth(args.input_prefix, genome_id, args.bam_folder, row[4], abs(int(row[3])-int(row[2])), args.samples.split(','))
                depth_per_pos[row[4]] = [overall, hp1, hp2, hp_all_soft, hp1_soft, hp2_soft]
    with open(genome_id+"/"+args.annotation_folder+"/"+genome_id+args.annotation_suffix, 'r') as in_an:
        for line in in_an:
            if count == 0:
                count = 1
                continue
            row = line.strip().split('\t')
            if count_per_pos[row[4]] < 2:
                truth_set["Novel"][row[1]] = row
            else:
                truth_set["HotSpot"][row[1]] = row
    #print(cut_start)
    #print(cut_end)
    # Check Pipeline first
#    control_inserts_per_chrom = defaultdict(IntervalTree)
#    with open(args.input_prefix+"/"+genome_id+"/Sample1/"+args.folder+"/Sample1.all.read_insertions.tsv", 'r') as in_data:
#        count = 0
#        for line in in_data:
#            if count == 0:
#                count += 1
#                continue
#            row = line.strip().split('\t')
#            control_inserts_per_chrom[row[0]][int(row[1])-500:int(row[2])+500] = 1
    for sample in args.samples.split(','):
        with open(args.input_prefix+"/"+genome_id+"/"+sample+"/"+args.folder+"/"+sample+args.suffix, 'r') as in_data:
            count = 0
            for line in in_data:
                if count == 0:
                    count += 1
                    continue
                row = line.strip().split('\t')
                #insert_list = control_inserts_per_chrom[row[0]][int(row[1]):int(row[2])]
                #if len(insert_list) == 0:
                if row[3] in truth_set[sample]:
                    #reads_to_check[row[3]] = 1
                    # have a read that we care about
                    #insert_size = truth_set[sample][row[3]][5]
                    #repeat_type = truth_set[sample][row[3]][6]
                    if row[8] == "PASS":
                        # Check more here
                        tp[row[3]] = 1#insert_size+"\t"+repeat_type+"\t"+str(read_depth[row[3]])
                    else:
                        #print(row)
                        hp = int(float(row[18]))
                        #print(hp)
                        depth_to_use = depth_per_pos[truth_set[sample][row[3]][4]][hp]
                        soft_count = depth_per_pos[truth_set[sample][row[3]][4]][hp+3]
                        if (count_per_pos[truth_set[sample][row[3]][4]]+soft_count)/depth_to_use > 0.3:
                            tn[row[3]] = 1
                        else:
                            fn[row[3]] = 1#insert_size+"\t"+repeat_type+"\t"+str(read_depth[row[3]])
                            #print(line.strip())
                            fn_rows[row[3]].append(row)
                else:
                    # Read we don't care about
                    if row[8] == "PASS":
                        # Check to see if this insert is actually an FP or if it occurs at a position where there is an insertion in the assembly
                        # BUT a position where we did not retain an insert supporting read for the experiement
                        # Mostly due to either the assembly insertion not being long enough, not mapping to repbase
                        # or the assembly contig having multiple mappings as it contained centromere sequence
                        if len(all_inserts[row[0]][int(row[1]):int(row[2])]) > 0:
                            continue
                        #insert_size = str(int(row[5])-int(row[4]))
                        #repeat_type = ""
                        #if "LINE" in line:
                        #    repeat_type = "LINE"
                        #elif "SINE" in line:
                        #    repeat_type = "SINE"
                        #elif "ERV" in line:
                        #    repeat_type = "ERV"
                        #elif "SVA" in line:
                        #    repeat_type = "SVA"
                        fp[row[3]] = 1
                        #print("\t".join(row))
                        #insert_size+"\t"+repeat_type
                        #reads_to_check[row[3]] = 1
                    else:
                        pass
                        #tn += 1
                #else:
                #    if row[3] in truth_set[sample]:
                #        insert_size = truth_set[sample][row[3]][5]
                #        repeat_type = truth_set[sample][row[3]][6]
                #        fn[row[3]] = insert_size+"\t"+repeat_type+"\t"+str(read_depth[row[3]])
                        #print(line.strip())
    #quit()
    # now check to see if we've missed anything.
    # Look for reads in truth set not found in
    true_fn = defaultdict(list)
    missed_fn = defaultdict(str)
    for sample in args.samples.split(','):
        for read in truth_set[sample]:
            if read in tp:
                #print(read)
                continue
            if read not in fn and read not in tn:
                #if (count_per_pos[truth_set[sample][read][4]]+depth_per_pos[truth_set[sample][read][4]][3])/depth_per_pos[truth_set[sample][read][4]][0] > 0.3:
                #    continue
                #insert_size = truth_set[sample][read][5]
                #repeat_type = truth_set[sample][read][6]
                missed_fn[read] = 1#insert_size+"\t"+repeat_type+"\t"+str(read_depth[read])
            #reads_to_check[read] = 1
            if read in fn:
                true_fn[read] = fn_rows[read]

    if args.fn:
        with open(args.fn, 'w') as out_fn:
            for read in true_fn:
                for row in true_fn[read]:
                    out_fn.write('\t'.join(row)+"\n")
        quit()
    # Remove any FN calls if we already have seen a TP for that read
    # Possible the read has more than one insertion and we have found the one we care abot
    # Other calls could be polymorphic or not RT inserts being called as FNs
    missed_fn_overall += len(missed_fn)
    combined_fn = defaultdict(str)
    for read in fn:
        if read not in tp:
            combined_fn[read] = fn[read]
    for read in missed_fn:
        combined_fn[read] = missed_fn[read]
    # Get the read_lengths
    #read_lengths = defaultdict(int)
    #for sample in ["Sample3", "Sample5"]:
    #    with pysam.FastaFile(args.input_prefix+"/"+genome_id+"/"+sample+"/fastq/"+sample+".fastq.gz") as fin:
    #        for read in reads_to_check:
    #            try:
    #                seq = fin.fetch(read)
    #                read_lengths[read] = len(seq)
    #            except KeyError:
    #                pass
    # Now figure out combined stats
    tp_novel = 0
    tp_hotspot = 0
    fn_novel = 0
    fn_hotspot = 0
    fn_novel_missed = 0
    fn_hotspot_missed = 0
    for read in tp:
        tp_overall += 1
        if read in truth_set["Novel"]:
            tp_novel += 1
        else:
            tp_hotspot +=1
        #insert_size = int(tp[read].split('\t')[0])
        #if insert_size < 500:
        #    tp_insert_size["Small"] += 1
        #elif insert_size < 2500:
        #    tp_insert_size["Medium"] += 1
        #elif insert_size < 10000:
        #    tp_insert_size["Large"] += 1
        #else:
        #    tp_insert_size["XLarge"] += 1
        #repeat_type = tp[read].split('\t')[1]
        #tp_repeat_type[repeat_type] += 1
        #if read_lengths[read] < 5000:
        #    tp_read_size["Small"] += 1
        #elif read_lengths[read] < 10000:
        #    tp_read_size["Medium"] += 1
        #elif read_lengths[read] < 15000:
        #    tp_read_size["Large"] += 1
        #else:
        #    tp_read_size["XLarge"] += 1
        #splice_depth = int(tp[read].split('\t')[2])
        #if splice_depth == 1:
        #    tp_spliced_depth["Single"] += 1
        #elif splice_depth == 2:
        #    tp_spliced_depth["Double"] += 1
        #elif splice_depth == 3:
        #    tp_spliced_depth["Triple"] += 1
        #elif splice_depth <= 5:
        #    tp_spliced_depth["Small"] += 1
        #elif splice_depth <= 10:
        #    tp_spliced_depth["Medium"] += 1
        #elif splice_depth <= 15:
        #    tp_spliced_depth["Large"] += 1
        #else:
        #    tp_spliced_depth["XLarge"] += 1
    for read in fp:
        fp_overall += 1
        #insert_size = int(fp[read].split('\t')[0])
        #if insert_size < 500:
        #    fp_insert_size["Small"] += 1
        #elif insert_size < 2500:
        #    fp_insert_size["Medium"] += 1
        #elif insert_size < 10000:
        #    fp_insert_size["Large"] += 1
        #else:
        #    fp_insert_size["XLarge"] += 1
        #repeat_type = fp[read].split('\t')[1]
        #fp_repeat_type[repeat_type] += 1
        #if read_lengths[read] < 5000:
        #    fp_read_size["Small"] += 1
        #elif read_lengths[read] < 10000:
        #    fp_read_size["Medium"] += 1
        #elif read_lengths[read] < 15000:
        #    fp_read_size["Large"] += 1
        #else:
        #    fp_read_size["XLarge"] += 1
    for read in combined_fn:
        fn_overall += 1
        if read in truth_set["Novel"]:
            fn_novel += 1
            if read in missed_fn:
                fn_novel_missed += 1
        else:
            fn_hotspot += 1
            if read in missed_fn:
                fn_hotspot_missed += 1
        #insert_size = int(combined_fn[read].split('\t')[0])
        #if insert_size < 500:
        #    fn_insert_size["Small"] += 1
        #elif insert_size < 2500:
        #    fn_insert_size["Medium"] += 1
        #elif insert_size < 10000:
        #    fn_insert_size["Large"] += 1
        #else:
        #    fn_insert_size["XLarge"] += 1
        #repeat_type = combined_fn[read].split('\t')[1]
        #fn_repeat_type[repeat_type] += 1
        #if read_lengths[read] < 5000:
        #    fn_read_size["Small"] += 1
        #elif read_lengths[read] < 10000:
        #    fn_read_size["Medium"] += 1
        #elif read_lengths[read] < 15000:
        #    fn_read_size["Large"] += 1
        #else:
        #    fn_read_size["XLarge"] += 1
        #splice_depth = int(combined_fn[read].split('\t')[2])
        #if splice_depth == 1:
        #    fn_spliced_depth["Single"] += 1
        #elif splice_depth == 2:
        #    fn_spliced_depth["Double"] += 1
        #elif splice_depth == 3:
        #    fn_spliced_depth["Triple"] += 1
        #elif splice_depth <= 5:
        #    fn_spliced_depth["Small"] += 1
        #elif splice_depth <= 10:
        #    fn_spliced_depth["Medium"] += 1
        #elif splice_depth <= 15:
        #    fn_spliced_depth["Large"] += 1
        #else:
        #    fn_spliced_depth["XLarge"] += 1

print("Pipeline:")
#tp_insert = []
#tp_reads = []
#tp_repeat = []
#tp_splice = []
#fp_insert = []
#fp_reads = []
#fp_repeat = []
#fn_insert = []
#fn_reads = []
#fn_repeat = []
#fn_splice = []
#for item in ["Small", "Medium", "Large", "XLarge"]:
#    tp_insert.append(str(tp_insert_size[item]))
#    tp_reads.append(str(tp_read_size[item]))
#    fp_insert.append(str(fp_insert_size[item]))
#    fp_reads.append(str(fp_read_size[item]))
#    fn_insert.append(str(fn_insert_size[item]))
#    fn_reads.append(str(fn_read_size[item]))
#for item in ["Single", "Double", "Triple", "Small", "Medium", "Large", "XLarge"]:
#    fn_splice.append(str(fn_spliced_depth[item]))
#    tp_splice.append(str(tp_spliced_depth[item]))
#for item in ["LINE", "SINE", "SVA", "ERV"]:
#    tp_repeat.append(str(tp_repeat_type[item]))
#    fp_repeat.append(str(fp_repeat_type[item]))
#    fn_repeat.append(str(fn_repeat_type[item]))
print("Overall: ")
print("TP : \t"+str(tp_overall))
print("FP : \t"+str(fp_overall))
print("FN : \t"+str(fn_overall))
print("Missed FN: \t"+str(missed_fn_overall))
print("TN: \t"+str(len(tn)))
print("Novel : \t"+str(tp_novel)+"\t"+str(fn_novel)+"\t"+str(fn_novel_missed))
print("HotSpot : \t"+str(tp_hotspot)+"\t"+str(fn_hotspot)+"\t"+str(fn_hotspot_missed))
#print("Insert Size: ")
#print("TP : \t"+"\t".join(tp_insert))
#print("FP : \t"+"\t".join(fp_insert))
#print("FN : \t"+"\t".join(fn_insert))
#print("Read Size: ")
#print("TP : \t"+"\t".join(tp_reads))
#print("FP : \t"+"\t".join(fp_reads))
#print("FN : \t"+"\t".join(fn_reads))
#print("Repeat Type: ")
#print("TP : \t"+"\t".join(tp_repeat))
#print("FP : \t"+"\t".join(fp_repeat))
#print("FN : \t"+"\t".join(fn_repeat))
#print("Spliced Depth: ")
#print("TP : \t"+"\t".join(tp_splice))
#print("FN : \t"+"\t".join(fn_splice))


#    # Now check Sniffles
#
#tp_overall = 0
#tp_insert_size = defaultdict(int)
#tp_read_size = defaultdict(int)
#tp_repeat_type = defaultdict(int)
#tp_spliced_depth = defaultdict(int)
#fp_insert_size = defaultdict(int)
#fp_overall = 0
#fp_read_size = defaultdict(int)
#fp_repeat_type = defaultdict(int)
#fp_spliced_depth = defaultdict(int)
#fn_insert_size = defaultdict(int)
#fn_read_size = defaultdict(int)
#fn_repeat_type = defaultdict(int)
#fn_spliced_depth = defaultdict(int)
#fn_overall = 0
#tn = 0
#
#tp_position_overall = 0
#tp_position_insert_size = defaultdict(int)
#tp_position_repeat_type = defaultdict(int)
#tp_position_spliced_depth = defaultdict(int)
#fp_position_insert_size = defaultdict(int)
#fp_position_overall = 0
#fp_position_repeat_type = defaultdict(int)
#fp_position_spliced_depth = defaultdict(int)
#fn_position_insert_size = defaultdict(int)
#fn_position_repeat_type = defaultdict(int)
#fn_position_spliced_depth = defaultdict(int)
#fn_position_overall = 0
#tn_position = 0
#
#for genome_id in args.id_list.split(','):
#    #print(genome_id)
#    reads_to_check = defaultdict(int)
#    tp = defaultdict(str)
#    fp = defaultdict(str)
#    fn = defaultdict(str)
#    tp_position = defaultdict(str)
#    fp_position = defaultdict(str)
#    fn_position = defaultdict(str)
#    # Open up the annotation file
#    truth_set = defaultdict(str)
#    truth_positions = defaultdict(IntervalTree)
#    reads_to_pos = defaultdict(list)
#    read_depth = defaultdict(int)
#    count = 0
#    with open(genome_id+"/"+args.annotation_folder+"/"+genome_id+args.annotation_suffix, 'r') as in_an:
#        for line in in_an:
#            if count == 0:
#                count = 1
#                continue
#            row = line.strip().split('\t')
#            truth_set[row[0]] = row[3]+"_"+row[4]
#            chrom = row[2].split(':')[0]
#            start = int(row[2].split(':')[1].split('-')[0])
#            end = int(row[2].split(':')[1].split('-')[1])
#            truth_positions[chrom][(start-500):(end+500)] = row[2]+"_"+row[3]+"_"+row[4]
#            reads_to_pos[row[2]].append(row[0])
#    for pos in reads_to_pos:
#        for read in reads_to_pos[pos]:
#            read_depth[read] = len(reads_to_pos[pos])
#    control_inserts_per_chrom = defaultdict(IntervalTree)
#    sniffles_control = pysam.VariantFile(args.input_prefix+"/"+genome_id+"/Sample1/"+args.sv_folder+"/Sample1"+args.sniffles_suffix)
#    for rec in sniffles_control.fetch():
#        if rec.info["SVTYPE"] == "INS" and int(rec.info["SVLEN"]) >= 50:
#            # Add this to the dict
#            control_inserts_per_chrom[rec.chrom][int(rec.start-500):int(rec.stop+500)] = 1
#    for sample in ["Sample3", "Sample5"]:
#        sniffles_sample = pysam.VariantFile(args.input_prefix+"/"+genome_id+"/"+sample+"/"+args.sv_folder+"/"+sample+args.sniffles_filtered_suffix)
#        for rec in sniffles_sample.fetch():
#            rec_row = str(rec).strip().split('\t')
#            if int(rec.info["SVLEN"]) >= 100:
#                insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
#                if len(insert_list) == 0 and rec_row[6] == "PASS":
#                    # Insert does not appear in Sample 1
#                    for name in rec.info["RNAMES"]:
#                        if name in truth_set:
#                            # Have a true positive here
#                            reads_to_check[name] = 1
#                            # have a read that we care about
#                            insert_size = truth_set[name].split('_')[0]
#                            repeat_type = truth_set[name].split('_')[1]
#                            tp[name] = insert_size+"\t"+repeat_type+"\t"+str(read_depth[name])
#                        else:
#                            repeat_type = rec.info["RTTYPE"]
#                            insert_size = str(rec.info["SVLEN"])
#                            fp[name] = insert_size+"\t"+repeat_type
#                            #reads_to_check[name] = 1
#                else:
#                    # Insertion expected but found to be polymorphic in Sample 1
#                    for name in rec.info["RNAMES"]:
#                        if name in truth_set:
#                            reads_to_check[name] = 1
#                            insert_size = truth_set[name].split('_')[0]
#                            repeat_type = truth_set[name].split('_')[1]
#                            fn[name] = insert_size+"\t"+repeat_type+"\t"+str(read_depth[name])
#                        else:
#                            tn += 1
#                insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
#                if len(insert_list) == 0 and rec_row[6] == "PASS":
#                    # Insert does not appear in Sample 1
#                    nearby = truth_positions[rec.chrom][rec.start:rec.stop]
#                    if len(nearby) > 0:
#                        for item in nearby:
#                            pos = item.data.split('_')[0]
#                            insert_size = item.data.split('_')[1]
#                            repeat_type = item.data.split('_')[2]
#                            tp_position[rec.chrom+":"+str(item.begin)+"-"+str(item.end)] = insert_size+"\t"+repeat_type+"\t"+str(len(reads_to_pos[pos]))
#                    else:
#                        repeat_type = rec.info["RTTYPE"]
#                        insert_size = str(rec.info["SVLEN"])
#                        fp_position[rec.chrom+":"+str(rec.start)+"-"+str(rec.stop)] = insert_size+"\t"+repeat_type
#                else:
#                    # Insertion expected but found to be polymorphic in Sample 1
#                    nearby = truth_positions[rec.chrom][rec.start:rec.stop]
#                    if len(nearby) > 0:
#                        for item in nearby:
#                            pos = item.data.split('_')[0]
#                            insert_size = item.data.split('_')[1]
#                            repeat_type = item.data.split('_')[2]
#                            fn_position[rec.chrom+":"+str(item.begin)+"-"+str(item.end)] = insert_size+"\t"+repeat_type+"\t"+str(len(reads_to_pos[pos]))
#                    else:
#                        tn_position += 1
#    # now check to see if we've missed anything.
#    # Look for reads in truth set not found in 
#    missed_fn = defaultdict(str)
#    missed_fn_position = defaultdict(str)
#    for read in truth_set:
#        if read in tp:
#            continue
#        if read not in fn:
#            insert_size = truth_set[read].split('_')[0]
#            repeat_type = truth_set[read].split('_')[1]
#            missed_fn[read] = insert_size+"\t"+repeat_type+"\t"+str(read_depth[read])
#            reads_to_check[read] = 1
#    for chrom in truth_positions:
#        for item in truth_positions[chrom]:
#            if chrom+":"+str(item.begin)+"-"+str(item.end) in tp_position:
#                continue
#            if chrom+":"+str(item.begin)+"-"+str(item.end) not in fn_position:
#                pos = item.data.split('_')[0]
#                insert_size = item.data.split('_')[1]
#                repeat_type = item.data.split('_')[2]
#                missed_fn_position[chrom+":"+str(item.begin)+"-"+str(item.end)] = insert_size+"\t"+repeat_type+"\t"+str(len(reads_to_pos[pos]))
#    # Remove any FN calls if we already have seen a TP for that read
#    # Possible the read has more than one insertion and we have found the one we care abot
#    # Other calls could be polymorphic or not RT inserts being called as FNs
#    combined_fn = defaultdict(str)
#    combined_fn_position = defaultdict(str)
#    for read in fn:
#        if read not in tp:
#            combined_fn[read] = fn[read]
#    for pos in fn_position:
#        if pos not in tp_position:
#            combined_fn_position[pos] = fn_position[pos]
#    for read in missed_fn:
#        combined_fn[read] = missed_fn[read]
#    for pos in missed_fn_position:
#        combined_fn_position[pos] = missed_fn_position[pos]
#    # Get the read_lengths
#    read_lengths = defaultdict(int)
#    #print(len(reads_to_check))
#    for sample in ["Sample3", "Sample5"]:
#        with pysam.FastaFile(args.input_prefix+"/"+genome_id+"/"+sample+"/fastq/"+sample+".fastq.gz") as fin:
#            for read in reads_to_check:
#                try:
#                    seq = fin.fetch(read)
#                    read_lengths[read] = len(seq)
#                except KeyError:
#                    pass
#    # Now figure out combined stats
#    for read in tp:
#        tp_overall += 1
#        insert_size = int(tp[read].split('\t')[0])
#        if insert_size < 500:
#            tp_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            tp_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            tp_insert_size["Large"] += 1
#        else:
#            tp_insert_size["XLarge"] += 1
#        repeat_type = tp[read].split('\t')[1]
#        tp_repeat_type[repeat_type] += 1
#        if read_lengths[read] < 5000:
#            tp_read_size["Small"] += 1
#        elif read_lengths[read] < 10000:
#            tp_read_size["Medium"] += 1
#        elif read_lengths[read] < 15000:
#            tp_read_size["Large"] += 1
#        else:
#            tp_read_size["XLarge"] += 1
#        splice_depth = int(tp[read].split('\t')[2])
#        if splice_depth == 1:
#            tp_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            tp_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            tp_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            tp_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            tp_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            tp_spliced_depth["Large"] += 1
#        else:
#            tp_spliced_depth["XLarge"] += 1
#    for pos in tp_position:
#        tp_position_overall += 1
#        insert_size = int(tp_position[pos].split('\t')[0])
#        if insert_size < 500:
#            tp_position_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            tp_position_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            tp_position_insert_size["Large"] += 1
#        else:
#            tp_position_insert_size["XLarge"] += 1
#        repeat_type = tp_position[pos].split('\t')[1]
#        tp_position_repeat_type[repeat_type] += 1
#        splice_depth = int(tp_position[pos].split('\t')[2])
#        if splice_depth == 1:
#            tp_position_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            tp_position_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            tp_position_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            tp_position_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            tp_position_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            tp_position_spliced_depth["Large"] += 1
#        else:
#            tp_position_spliced_depth["XLarge"] += 1
#    for read in fp:
#        fp_overall += 1
#        insert_size = int(fp[read].split('\t')[0])
#        if insert_size < 500:
#            fp_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fp_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fp_insert_size["Large"] += 1
#        else:
#            fp_insert_size["XLarge"] += 1
#        repeat_type = fp[read].split('\t')[1]
#        fp_repeat_type[repeat_type] += 1
#        if read_lengths[read] < 5000:
#            fp_read_size["Small"] += 1
#        elif read_lengths[read] < 10000:
#            fp_read_size["Medium"] += 1
#        elif read_lengths[read] < 15000:
#            fp_read_size["Large"] += 1
#        else:
#            fp_read_size["XLarge"] += 1
#    for pos in fp_position:
#        fp_position_overall += 1
#        insert_size = int(fp_position[pos].split('\t')[0])
#        if insert_size < 500:
#            fp_position_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fp_position_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fp_position_insert_size["Large"] += 1
#        else:
#            fp_position_insert_size["XLarge"] += 1
#        repeat_type = fp_position[pos].split('\t')[1]
#        fp_position_repeat_type[repeat_type] += 1
#    for read in combined_fn:
#        fn_overall += 1
#        insert_size = int(combined_fn[read].split('\t')[0])
#        if insert_size < 500:
#            fn_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fn_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fn_insert_size["Large"] += 1
#        else:
#            fn_insert_size["XLarge"] += 1
#        repeat_type = combined_fn[read].split('\t')[1]
#        fn_repeat_type[repeat_type] += 1
#        if read_lengths[read] < 5000:
#            fn_read_size["Small"] += 1
#        elif read_lengths[read] < 10000:
#            fn_read_size["Medium"] += 1
#        elif read_lengths[read] < 15000:
#            fn_read_size["Large"] += 1
#        else:
#            fn_read_size["XLarge"] += 1
#        splice_depth = int(combined_fn[read].split('\t')[2])
#        if splice_depth == 1:
#            fn_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            fn_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            fn_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            fn_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            fn_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            fn_spliced_depth["Large"] += 1
#        else:
#            fn_spliced_depth["XLarge"] += 1
#    for pos in combined_fn_position:
#        fn_position_overall += 1
#        insert_size = int(combined_fn_position[pos].split('\t')[0])
#        if insert_size < 500:
#            fn_position_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fn_position_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fn_position_insert_size["Large"] += 1
#        else:
#            fn_position_insert_size["XLarge"] += 1
#        repeat_type = combined_fn_position[pos].split('\t')[1]
#        fn_position_repeat_type[repeat_type] += 1
#        splice_depth = int(combined_fn_position[pos].split('\t')[2])
#        if splice_depth == 1:
#            fn_position_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            fn_position_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            fn_position_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            fn_position_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            fn_position_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            fn_position_spliced_depth["Large"] += 1
#        else:
#            fn_position_spliced_depth["XLarge"] += 1
#
#print("Sniffles:")
#tp_insert = []
#tp_reads = []
#tp_repeat = []
#tp_splice = []
#fp_insert = []
#fp_reads = []
#fp_repeat = []
#fn_insert = []
#fn_reads = []
#fn_repeat = []
#fn_splice = []
#for item in ["Small", "Medium", "Large", "XLarge"]:
#    tp_insert.append(str(tp_insert_size[item]))
#    tp_reads.append(str(tp_read_size[item]))
#    fp_insert.append(str(fp_insert_size[item]))
#    fp_reads.append(str(fp_read_size[item]))
#    fn_insert.append(str(fn_insert_size[item]))
#    fn_reads.append(str(fn_read_size[item]))
#for item in ["Single", "Double", "Triple", "Small", "Medium", "Large", "XLarge"]:
#    fn_splice.append(str(fn_spliced_depth[item]))
#    tp_splice.append(str(tp_spliced_depth[item]))
#for item in ["LINE", "SINE", "SVA", "ERV"]:
#    tp_repeat.append(str(tp_repeat_type[item]))
#    fp_repeat.append(str(fp_repeat_type[item]))
#    fn_repeat.append(str(fn_repeat_type[item]))
#print("Overall: ")
#print("TP : \t"+str(tp_overall))
#print("FP : \t"+str(fp_overall))
#print("FN : \t"+str(fn_overall))
#print("Insert Size: ")
#print("TP : \t"+"\t".join(tp_insert))
#print("FP : \t"+"\t".join(fp_insert))
#print("FN : \t"+"\t".join(fn_insert))
#print("Read Size: ")
#print("TP : \t"+"\t".join(tp_reads))
#print("FP : \t"+"\t".join(fp_reads))
#print("FN : \t"+"\t".join(fn_reads))
#print("Repeat Type: ")
#print("TP : \t"+"\t".join(tp_repeat))
#print("FP : \t"+"\t".join(fp_repeat))
#print("FN : \t"+"\t".join(fn_repeat))
#print("Spliced Depth: ")
#print("TP : \t"+"\t".join(tp_splice))
#print("FN : \t"+"\t".join(fn_splice))
#
#print("Sniffles Positions:")
#tp_insert = []
#tp_repeat = []
#tp_splice = []
#fp_insert = []
#fp_repeat = []
#fn_insert = []
#fn_repeat = []
#fn_splice = []
#for item in ["Small", "Medium", "Large", "XLarge"]:
#    tp_insert.append(str(tp_position_insert_size[item]))
#    fp_insert.append(str(fp_position_insert_size[item]))
#    fn_insert.append(str(fn_position_insert_size[item]))
#for item in ["Single", "Double", "Triple", "Small", "Medium", "Large", "XLarge"]:
#    fn_splice.append(str(fn_position_spliced_depth[item]))
#    tp_splice.append(str(tp_position_spliced_depth[item]))
#for item in ["LINE", "SINE", "SVA", "ERV"]:
#    tp_repeat.append(str(tp_position_repeat_type[item]))
#    fp_repeat.append(str(fp_position_repeat_type[item]))
#    fn_repeat.append(str(fn_position_repeat_type[item]))
#print("Overall: ")
#print("TP : \t"+str(tp_position_overall))
#print("FP : \t"+str(fp_position_overall))
#print("FN : \t"+str(fn_position_overall))
#print("Insert Size: ")
#print("TP : \t"+"\t".join(tp_insert))
#print("FP : \t"+"\t".join(fp_insert))
#print("FN : \t"+"\t".join(fn_insert))
#print("Repeat Type: ")
#print("TP : \t"+"\t".join(tp_repeat))
#print("FP : \t"+"\t".join(fp_repeat))
#print("FN : \t"+"\t".join(fn_repeat))
#print("Spliced Depth: ")
#print("TP : \t"+"\t".join(tp_splice))
#print("FN : \t"+"\t".join(fn_splice))
#
#
##    # Now check CuteSV
#
#tp_overall = 0
#tp_insert_size = defaultdict(int)
#tp_read_size = defaultdict(int)
#tp_repeat_type = defaultdict(int)
#tp_spliced_depth = defaultdict(int)
#fp_insert_size = defaultdict(int)
#fp_overall = 0
#fp_read_size = defaultdict(int)
#fp_repeat_type = defaultdict(int)
#fp_spliced_depth = defaultdict(int)
#fn_insert_size = defaultdict(int)
#fn_read_size = defaultdict(int)
#fn_repeat_type = defaultdict(int)
#fn_spliced_depth = defaultdict(int)
#fn_overall = 0
#tn = 0
#
#tp_position_overall = 0
#tp_position_insert_size = defaultdict(int)
#tp_position_repeat_type = defaultdict(int)
#tp_position_spliced_depth = defaultdict(int)
#fp_position_insert_size = defaultdict(int)
#fp_position_overall = 0
#fp_position_repeat_type = defaultdict(int)
#fp_position_spliced_depth = defaultdict(int)
#fn_position_insert_size = defaultdict(int)
#fn_position_repeat_type = defaultdict(int)
#fn_position_spliced_depth = defaultdict(int)
#fn_position_overall = 0
#tn_position = 0
#
#for genome_id in args.id_list.split(','):
#    #print(genome_id)
#    reads_to_check = defaultdict(int)
#    tp = defaultdict(str)
#    fp = defaultdict(str)
#    fn = defaultdict(str)
#    tp_position = defaultdict(str)
#    fp_position = defaultdict(str)
#    fn_position = defaultdict(str)
#    # Open up the annotation file
#    truth_set = defaultdict(str)
#    truth_positions = defaultdict(IntervalTree)
#    reads_to_pos = defaultdict(list)
#    read_depth = defaultdict(int)
#    count = 0
#    with open(genome_id+"/"+args.annotation_folder+"/"+genome_id+args.annotation_suffix, 'r') as in_an:
#        for line in in_an:
#            if count == 0:
#                count = 1
#                continue
#            row = line.strip().split('\t')
#            truth_set[row[0]] = row[3]+"_"+row[4]
#            chrom = row[2].split(':')[0]
#            start = int(row[2].split(':')[1].split('-')[0])
#            end = int(row[2].split(':')[1].split('-')[1])
#            truth_positions[chrom][(start-500):(end+500)] = row[2]+"_"+row[3]+"_"+row[4]
#            reads_to_pos[row[2]].append(row[0])
#    print(len(truth_set))
#    for pos in reads_to_pos:
#        for read in reads_to_pos[pos]:
#            read_depth[read] = len(reads_to_pos[pos])
#    control_inserts_per_chrom = defaultdict(IntervalTree)
#    sniffles_control = pysam.VariantFile(args.input_prefix+"/"+genome_id+"/Sample1/"+args.sv_folder+"/Sample1"+args.cutesv_suffix)
#    for rec in sniffles_control.fetch():
#        if rec.info["SVTYPE"] == "INS" and int(rec.info["SVLEN"]) >= 50:
#            # Add this to the dict
#            control_inserts_per_chrom[rec.chrom][int(rec.start-500):int(rec.stop+500)] = 1
#    for sample in ["Sample3", "Sample5"]:
#        sniffles_sample = pysam.VariantFile(args.input_prefix+"/"+genome_id+"/"+sample+"/"+args.sv_folder+"/"+sample+args.cutesv_filtered_suffix)
#        for rec in sniffles_sample.fetch():
#            rec_row = str(rec).strip().split('\t')
#            if int(rec.info["SVLEN"]) >= 100:
#                insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
#                if len(insert_list) == 0 and rec_row[6] == "PASS":
#                    # Insert does not appear in Sample 1
#                    for name in rec.info["RNAMES"]:
#                        if name in truth_set:
#                            # Have a true positive here
#                            reads_to_check[name] = 1
#                            # have a read that we care about
#                            insert_size = truth_set[name].split('_')[0]
#                            repeat_type = truth_set[name].split('_')[1]
#                            tp[name] = insert_size+"\t"+repeat_type+"\t"+str(read_depth[name])
#                        else:
#                            repeat_type = rec.info["RTTYPE"]
#                            insert_size = str(rec.info["SVLEN"])
#                            fp[name] = insert_size+"\t"+repeat_type
#                            #reads_to_check[name] = 1
#                else:
#                    # Insertion expected but found to be polymorphic in Sample 1
#                    for name in rec.info["RNAMES"]:
#                        if name in truth_set:
#                            reads_to_check[name] = 1
#                            insert_size = truth_set[name].split('_')[0]
#                            repeat_type = truth_set[name].split('_')[1]
#                            fn[name] = insert_size+"\t"+repeat_type+"\t"+str(read_depth[name])
#                        else:
#                            tn += 1
#                insert_list = control_inserts_per_chrom[rec.chrom][rec.start:rec.stop]
#                if len(insert_list) == 0 and rec_row[6] == "PASS":
#                    # Insert does not appear in Sample 1
#                    nearby = truth_positions[rec.chrom][rec.start:rec.stop]
#                    if len(nearby) > 0:
#                        for item in nearby:
#                            pos = item.data.split('_')[0]
#                            insert_size = item.data.split('_')[1]
#                            repeat_type = item.data.split('_')[2]
#                            tp_position[rec.chrom+":"+str(item.begin)+"-"+str(item.end)] = insert_size+"\t"+repeat_type+"\t"+str(len(reads_to_pos[pos]))
#                    else:
#                        repeat_type = rec.info["RTTYPE"]
#                        insert_size = str(rec.info["SVLEN"])
#                        fp_position[rec.chrom+":"+str(rec.start)+"-"+str(rec.stop)] = insert_size+"\t"+repeat_type
#                else:
#                    # Insertion expected but found to be polymorphic in Sample 1
#                    nearby = truth_positions[rec.chrom][rec.start:rec.stop]
#                    if len(nearby) > 0:
#                        for item in nearby:
#                            pos = item.data.split('_')[0]
#                            insert_size = item.data.split('_')[1]
#                            repeat_type = item.data.split('_')[2]
#                            fn_position[rec.chrom+":"+str(item.begin)+"-"+str(item.end)] = insert_size+"\t"+repeat_type+"\t"+str(len(reads_to_pos[pos]))
#                    else:
#                        tn_position += 1
#    # now check to see if we've missed anything.
#    # Look for reads in truth set not found in 
#    missed_fn = defaultdict(str)
#    missed_fn_position = defaultdict(str)
#    for read in truth_set:
#        if read in tp:
#            continue
#        if read not in fn:
#            insert_size = truth_set[read].split('_')[0]
#            repeat_type = truth_set[read].split('_')[1]
#            missed_fn[read] = insert_size+"\t"+repeat_type+"\t"+str(read_depth[read])
#            reads_to_check[read] = 1
#    for chrom in truth_positions:
#        for item in truth_positions[chrom]:
#            if chrom+":"+str(item.begin)+"-"+str(item.end) in tp_position:
#                continue
#            if chrom+":"+str(item.begin)+"-"+str(item.end) not in fn_position:
#                pos = item.data.split('_')[0]
#                insert_size = item.data.split('_')[1]
#                repeat_type = item.data.split('_')[2]
#                missed_fn_position[chrom+":"+str(item.begin)+"-"+str(item.end)] = insert_size+"\t"+repeat_type+"\t"+str(len(reads_to_pos[pos]))
#    # Remove any FN calls if we already have seen a TP for that read
#    # Possible the read has more than one insertion and we have found the one we care abot
#    # Other calls could be polymorphic or not RT inserts being called as FNs
#    combined_fn = defaultdict(str)
#    combined_fn_position = defaultdict(str)
#    for read in fn:
#        if read not in tp:
#            combined_fn[read] = fn[read]
#    for pos in fn_position:
#        if pos not in tp_position:
#            combined_fn_position[pos] = fn_position[pos]
#    for read in missed_fn:
#        combined_fn[read] = missed_fn[read]
#    for pos in missed_fn_position:
#        combined_fn_position[pos] = missed_fn_position[pos]
#    # Get the read_lengths
#    read_lengths = defaultdict(int)
#    #print(len(reads_to_check))
#    for sample in ["Sample3", "Sample5"]:
#        with pysam.FastaFile(args.input_prefix+"/"+genome_id+"/"+sample+"/fastq/"+sample+".fastq.gz") as fin:
#            for read in reads_to_check:
#                try:
#                    seq = fin.fetch(read)
#                    read_lengths[read] = len(seq)
#                except KeyError:
#                    pass
#    # Now figure out combined stats
#    for read in tp:
#        tp_overall += 1
#        insert_size = int(tp[read].split('\t')[0])
#        if insert_size < 500:
#            tp_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            tp_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            tp_insert_size["Large"] += 1
#        else:
#            tp_insert_size["XLarge"] += 1
#        repeat_type = tp[read].split('\t')[1]
#        tp_repeat_type[repeat_type] += 1
#        if read_lengths[read] < 5000:
#            tp_read_size["Small"] += 1
#        elif read_lengths[read] < 10000:
#            tp_read_size["Medium"] += 1
#        elif read_lengths[read] < 15000:
#            tp_read_size["Large"] += 1
#        else:
#            tp_read_size["XLarge"] += 1
#        splice_depth = int(tp[read].split('\t')[2])
#        if splice_depth == 1:
#            tp_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            tp_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            tp_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            tp_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            tp_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            tp_spliced_depth["Large"] += 1
#        else:
#            tp_spliced_depth["XLarge"] += 1
#    for pos in tp_position:
#        tp_position_overall += 1
#        insert_size = int(tp_position[pos].split('\t')[0])
#        if insert_size < 500:
#            tp_position_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            tp_position_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            tp_position_insert_size["Large"] += 1
#        else:
#            tp_position_insert_size["XLarge"] += 1
#        repeat_type = tp_position[pos].split('\t')[1]
#        tp_position_repeat_type[repeat_type] += 1
#        splice_depth = int(tp_position[pos].split('\t')[2])
#        if splice_depth == 1:
#            tp_position_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            tp_position_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            tp_position_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            tp_position_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            tp_position_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            tp_position_spliced_depth["Large"] += 1
#        else:
#            tp_position_spliced_depth["XLarge"] += 1
#    for read in fp:
#        fp_overall += 1
#        insert_size = int(fp[read].split('\t')[0])
#        if insert_size < 500:
#            fp_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fp_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fp_insert_size["Large"] += 1
#        else:
#            fp_insert_size["XLarge"] += 1
#        repeat_type = fp[read].split('\t')[1]
#        fp_repeat_type[repeat_type] += 1
#        if read_lengths[read] < 5000:
#            fp_read_size["Small"] += 1
#        elif read_lengths[read] < 10000:
#            fp_read_size["Medium"] += 1
#        elif read_lengths[read] < 15000:
#            fp_read_size["Large"] += 1
#        else:
#            fp_read_size["XLarge"] += 1
#    for pos in fp_position:
#        fp_position_overall += 1
#        insert_size = int(fp_position[pos].split('\t')[0])
#        if insert_size < 500:
#            fp_position_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fp_position_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fp_position_insert_size["Large"] += 1
#        else:
#            fp_position_insert_size["XLarge"] += 1
#        repeat_type = fp_position[pos].split('\t')[1]
#        fp_position_repeat_type[repeat_type] += 1
#    for read in combined_fn:
#        fn_overall += 1
#        insert_size = int(combined_fn[read].split('\t')[0])
#        if insert_size < 500:
#            fn_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fn_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fn_insert_size["Large"] += 1
#        else:
#            fn_insert_size["XLarge"] += 1
#        repeat_type = combined_fn[read].split('\t')[1]
#        fn_repeat_type[repeat_type] += 1
#        if read_lengths[read] < 5000:
#            fn_read_size["Small"] += 1
#        elif read_lengths[read] < 10000:
#            fn_read_size["Medium"] += 1
#        elif read_lengths[read] < 15000:
#            fn_read_size["Large"] += 1
#        else:
#            fn_read_size["XLarge"] += 1
#        splice_depth = int(combined_fn[read].split('\t')[2])
#        if splice_depth == 1:
#            fn_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            fn_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            fn_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            fn_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            fn_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            fn_spliced_depth["Large"] += 1
#        else:
#            fn_spliced_depth["XLarge"] += 1
#    for pos in combined_fn_position:
#        fn_position_overall += 1
#        insert_size = int(combined_fn_position[pos].split('\t')[0])
#        if insert_size < 500:
#            fn_position_insert_size["Small"] += 1
#        elif insert_size < 2500:
#            fn_position_insert_size["Medium"] += 1
#        elif insert_size < 10000:
#            fn_position_insert_size["Large"] += 1
#        else:
#            fn_position_insert_size["XLarge"] += 1
#        repeat_type = combined_fn_position[pos].split('\t')[1]
#        fn_position_repeat_type[repeat_type] += 1
#        splice_depth = int(combined_fn_position[pos].split('\t')[2])
#        if splice_depth == 1:
#            fn_position_spliced_depth["Single"] += 1
#        elif splice_depth == 2:
#            fn_position_spliced_depth["Double"] += 1
#        elif splice_depth == 3:
#            fn_position_spliced_depth["Triple"] += 1
#        elif splice_depth <= 5:
#            fn_position_spliced_depth["Small"] += 1
#        elif splice_depth <= 10:
#            fn_position_spliced_depth["Medium"] += 1
#        elif splice_depth <= 15:
#            fn_position_spliced_depth["Large"] += 1
#        else:
#            fn_position_spliced_depth["XLarge"] += 1
#
#print("CuteSV:")
#tp_insert = []
#tp_reads = []
#tp_repeat = []
#tp_splice = []
#fp_insert = []
#fp_reads = []
#fp_repeat = []
#fn_insert = []
#fn_reads = []
#fn_repeat = []
#fn_splice = []
#for item in ["Small", "Medium", "Large", "XLarge"]:
#    tp_insert.append(str(tp_insert_size[item]))
#    tp_reads.append(str(tp_read_size[item]))
#    fp_insert.append(str(fp_insert_size[item]))
#    fp_reads.append(str(fp_read_size[item]))
#    fn_insert.append(str(fn_insert_size[item]))
#    fn_reads.append(str(fn_read_size[item]))
#for item in ["Single", "Double", "Triple", "Small", "Medium", "Large", "XLarge"]:
#    fn_splice.append(str(fn_spliced_depth[item]))
#    tp_splice.append(str(tp_spliced_depth[item]))
#for item in ["LINE", "SINE", "SVA", "ERV"]:
#    tp_repeat.append(str(tp_repeat_type[item]))
#    fp_repeat.append(str(fp_repeat_type[item]))
#    fn_repeat.append(str(fn_repeat_type[item]))
#print("Overall: ")
#print("TP : \t"+str(tp_overall))
#print("FP : \t"+str(fp_overall))
#print("FN : \t"+str(fn_overall))
#print("Insert Size: ")
#print("TP : \t"+"\t".join(tp_insert))
#print("FP : \t"+"\t".join(fp_insert))
#print("FN : \t"+"\t".join(fn_insert))
#print("Read Size: ")
#print("TP : \t"+"\t".join(tp_reads))
#print("FP : \t"+"\t".join(fp_reads))
#print("FN : \t"+"\t".join(fn_reads))
#print("Repeat Type: ")
#print("TP : \t"+"\t".join(tp_repeat))
#print("FP : \t"+"\t".join(fp_repeat))
#print("FN : \t"+"\t".join(fn_repeat))
#print("Spliced Depth: ")
#print("TP : \t"+"\t".join(tp_splice))
#print("FN : \t"+"\t".join(fn_splice))
#
#print("CuteSV Positions:")
#tp_insert = []
#tp_repeat = []
#tp_splice = []
#fp_insert = []
#fp_repeat = []
#fn_insert = []
#fn_repeat = []
#fn_splice = []
#for item in ["Small", "Medium", "Large", "XLarge"]:
#    tp_insert.append(str(tp_position_insert_size[item]))
#    fp_insert.append(str(fp_position_insert_size[item]))
#    fn_insert.append(str(fn_position_insert_size[item]))
#for item in ["Single", "Double", "Triple", "Small", "Medium", "Large", "XLarge"]:
#    fn_splice.append(str(fn_position_spliced_depth[item]))
#    tp_splice.append(str(tp_position_spliced_depth[item]))
#for item in ["LINE", "SINE", "SVA", "ERV"]:
#    tp_repeat.append(str(tp_position_repeat_type[item]))
#    fp_repeat.append(str(fp_position_repeat_type[item]))
#    fn_repeat.append(str(fn_position_repeat_type[item]))
#print("Overall: ")
#print("TP : \t"+str(tp_position_overall))
#print("FP : \t"+str(fp_position_overall))
#print("FN : \t"+str(fn_position_overall))
#print("Insert Size: ")
#print("TP : \t"+"\t".join(tp_insert))
#print("FP : \t"+"\t".join(fp_insert))
#print("FN : \t"+"\t".join(fn_insert))
#print("Repeat Type: ")
#print("TP : \t"+"\t".join(tp_repeat))
#print("FP : \t"+"\t".join(fp_repeat))
#print("FN : \t"+"\t".join(fn_repeat))
#print("Spliced Depth: ")
#print("TP : \t"+"\t".join(tp_splice))
#print("FN : \t"+"\t".join(fn_splice))
#
#
