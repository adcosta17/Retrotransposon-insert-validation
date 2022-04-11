import os
import glob

rule all_assembly_bams:
    input:
        expand("{s}/assembly_mapped/{s}.sorted.bam.bai", s=config["samples"])

rule all_assembly_insertions:
	input:
		expand("{s}/assembly_analysis/{s}.insertions.tsv", s=config["samples"])

rule all_assembly_annotated_insertions:
	input:
		expand("{s}/assembly_analysis/{s}.insertions.repbase_annotated.tsv", s=config["samples"])

rule all_assembly_combined:
	input:
		"combined_assembly/combined_multi_sample.tsv"

rule all_read_assembly_mapped:
	input:
		expand("{s}/reads_assembly_mapped/{s}.sorted.1.paf.gz", s=config["samples"])

rule all_insert_supporting_reads:
	input:
		expand("{s}/reads_assembly_mapped/{s}.insert_supporting_reads.tsv", s=config["samples"])

rule all_modified_fastq:
	input:
		expand("{s}/fastq_inserted/{s}.insert_annotation.tsv", s=config["samples"])

rule all_metrics:
	input:
		expand("{s}/metrics/{s}.metrics.tsv", s=config["samples"])

rule all_short_read_fastq:
	input:
		expand("{s}/illumina/{s}.p1.fastq.gz", s=config["samples_subset"])

rule all_abyss_assemblies:
	input:
		expand("{s}/illumina_assembly/assembly-scaffolds.fa", s=config["samples_subset"])

rule all_downsampled_fastq:
	input:
		expand("{s}/illumina/{s}.p1.0.5.fastq.gz", s=config["samples_subset"]),
        expand("{s}/illumina/{s}.p2.0.5.fastq.gz", s=config["samples_subset"])

rule all_extracted_insertions:
	input:
		expand("{s}/reads_assembly_mapped/{s}.extracted_insertions.tsv",s=config["samples"])

rule all_insert_added_fastqs:
	input:
		expand("{s}/fastq_insert_added/{s}.insert_annotation.tsv",s=config["samples"])

rule all_insert_added_position_fastqs:
	input:
		expand("{s}/fastq_insert_added_position/{s}.insert_annotation.tsv",s=config["samples"])

rule all_mat_pat_bams:
	input:
		expand("{s}/assembly_mapped/{s}.mat.sorted.bam",s=config["samples"]),
		expand("{s}/assembly_mapped/{s}.pat.sorted.bam",s=config["samples"])

rule all_mat_pat_read_bams:
	input:
		expand("{s}/reads_assembly_mapped/{s}.mat.sorted.1.bam",s=config["samples"]),
        expand("{s}/reads_assembly_mapped/{s}.pat.sorted.1.bam",s=config["samples"]),
        expand("{s}/reads_assembly_mapped/{s}.mat.sorted.1.bam.bai",s=config["samples"]),
        expand("{s}/reads_assembly_mapped/{s}.pat.sorted.1.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_mapped/{s}.mat.sorted.5.bam",s=config["samples"]),
        expand("{s}/reads_assembly_mapped/{s}.mat.sorted.3.bam",s=config["samples"]),
        expand("{s}/reads_assembly_mapped/{s}.pat.sorted.5.bam",s=config["samples"]),
        expand("{s}/reads_assembly_mapped/{s}.pat.sorted.3.bam",s=config["samples"])

rule all_mat_pat_subset_to_ref_inserts:
	input:
		expand("{s}/assembly_subset/{s}.mat.insertions.repbase_annotated.tsv",s=config["samples"]),
        expand("{s}/assembly_subset/{s}.pat.insertions.repbase_annotated.tsv",s=config["samples"])

rule all_extracted_insert_positions:
	input:
		expand("{s}/reads_assembly_mapped/{s}.extracted_insertions_with_positions.tsv",s=config["samples"])

rule all_mat_mapped_to_pat_bams:
	input:
		expand("{s}/pat_to_mat/{s}.pat_to_mat.bam",s=config["samples"])

rule all_mat_mapped_to_pat_bams_plus_fastq:
	input:
		expand("{s}/pat_to_mat/{s}.pat_to_mat.bam",s=config["samples"]),
		expand("{s}/pat_to_mat/{s}.pat_to_mat.bam.bai",s=config["samples"]),
		expand("{s}/fastq_insert_added_position/{s}.insert_annotation.tsv",s=config["samples"])

rule all_mat_to_pat_inserts:
	input:
		expand("{s}/pat_to_mat/{s}.insertions.repbase_annotated.tsv",s=config["samples"])

rule all_reads_mat_inserts:
	input:
		expand("{s}/pat_to_mat/{s}.reads.1.insertions.tsv", s=config["samples"]),
		expand("{s}/pat_to_mat/{s}.reads.3.insertions.tsv", s=config["samples"]),
		expand("{s}/pat_to_mat/{s}.reads.5.insertions.tsv", s=config["samples"])

rule all_deletion_fastqs:
	input:
		expand("{s}/fastq_with_deletions/{s}.annotation.tsv",s=config["samples"]),
		expand("{s}/fastq_with_deletions/Sample1.fastq.gz",s=config["samples"]),
		expand("{s}/fastq_with_deletions/Sample3.fastq.gz",s=config["samples"]),
		expand("{s}/fastq_with_deletions/Sample5.fastq.gz",s=config["samples"])


rule all_reads_to_subset_bams:
	input:
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.1.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.3.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.5.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.1.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.3.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.5.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.1.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.3.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.5.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.1.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.3.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.5.bam.bai",s=config["samples"])

rule all_reads_to_del_bams:
	input:
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.1.deletions.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.3.deletions.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.5.deletions.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.1.deletions.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.3.deletions.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.5.deletions.bam",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.1.deletions.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.3.deletions.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.sorted.5.deletions.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.1.deletions.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.3.deletions.bam.bai",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.mat.sorted.5.deletions.bam.bai",s=config["samples"])

rule all_mtp_and_ptm:
	input:
		expand("{s}/pat_to_mat/{s}.pat_subset_to_mat_subset.bam",s=config["samples"]),
		expand("{s}/mat_to_pat/{s}.mat_subset_to_pat_subset.bam",s=config["samples"]),
		expand("{s}/pat_to_mat/{s}.pat_subset_to_mat_subset.bam.bai",s=config["samples"]),
		expand("{s}/mat_to_pat/{s}.mat_subset_to_pat_subset.bam.bai",s=config["samples"])

rule all_reads_mtp_ptm:
	input:
		expand("{s}/reads_assembly_subset_mapped/{s}.pat.filtered.sorted.1.bam",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.pat.filtered.sorted.1.bam.bai",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.pat.filtered.sorted.3.bam",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.pat.filtered.sorted.3.bam.bai",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.pat.filtered.sorted.5.bam",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.pat.filtered.sorted.5.bam.bai",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.filtered.sorted.1.bam",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.filtered.sorted.1.bam.bai",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.filtered.sorted.3.bam",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.filtered.sorted.3.bam.bai",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.filtered.sorted.5.bam",s=config["samples"]),
        expand("{s}/reads_assembly_subset_mapped/{s}.mat.filtered.sorted.5.bam.bai",s=config["samples"])

rule all_insert_reads_mtp_ptm:
	input:
		expand("{s}/reads_assembly_subset_mapped/{s}.reads.pat.1.insertions.tsv",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.reads.pat.3.insertions.tsv",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.reads.pat.5.insertions.tsv",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.reads.mat.1.insertions.tsv",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.reads.mat.3.insertions.tsv",s=config["samples"]),
		expand("{s}/reads_assembly_subset_mapped/{s}.reads.mat.5.insertions.tsv",s=config["samples"])

rule all_deletion_fastq_and_alignments:
	input:
		expand("{s}/fastq_with_deletions/{s}.annotation.tsv",s=config["samples"]),
		expand("{s}/fastq_with_deletions/Sample1.fastq.gz",s=config["samples"]),
		expand("{s}/fastq_with_deletions/Sample3.fastq.gz",s=config["samples"]),
		expand("{s}/fastq_with_deletions/Sample5.fastq.gz",s=config["samples"]),
		expand("{s}/assembly_subset/{s}.mat.insertions.repbase_annotated.tsv",s=config["samples"]),
        expand("{s}/assembly_subset/{s}.pat.insertions.repbase_annotated.tsv",s=config["samples"]),
        expand("{s}/assembly_analysis/{s}.insertions.repbase_annotated.tsv", s=config["samples"])


#rule all_mapped_other_assembly:
#	input:
#		expand("{s}/{s}.reads_mapped_to.{o}.mat.bam",s=config["samples_plus_extra"],o=config["samples"]),
#        expand("{s}/{s}.reads_mapped_to.{o}.pat.bam",s=config["samples_plus_extra"],o=config["samples"]),
#        expand("{s}/{s}.reads_mapped_to.grch38.bam",s=config["samples_plus_extra"]),
#        expand("{s}/{s}.reads_mapped_to.combined.bam",s=config["samples_plus_extra"]),
#        expand("{s}/{s}.reads_mapped_to.{o}.mat.bam.bai",s=config["samples_plus_extra"],o=config["samples"]),
#        expand("{s}/{s}.reads_mapped_to.{o}.pat.bam.bai",s=config["samples_plus_extra"],o=config["samples"]),
#        expand("{s}/{s}.reads_mapped_to.grch38.bam.bai",s=config["samples_plus_extra"]),
#        expand("{s}/{s}.reads_mapped_to.combined.bam.bai",s=config["samples_plus_extra"])


rule all_hg0002_deletion_fastq:
	input:
		expand("{s}/fastq_with_deletions/HG0002.{i}.fastq.gz", i=config["counts"], s=config["samples"])

include: "rules/mapping.smk"
include: "rules/rt_winnow_validation.smk"
include: "rules/short_read.smk"
include: "rules/deletion_validation.smk"
include: "rules/hg0002.smk"

