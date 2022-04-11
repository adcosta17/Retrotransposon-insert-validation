
# Validation rules

def get_read_first_fastq(wildcards):
    return wildcards.sample+"/fastq/"+wildcards.sample+"_1_Guppy_4.0.11_prom.fastq.gz"

def get_read_third_fastq(wildcards):
    return wildcards.sample+"/fastq/"+wildcards.sample+"_3_Guppy_4.0.11_prom.fastq.gz"

def get_read_fifth_fastq(wildcards):
    return wildcards.sample+"/fastq/"+wildcards.sample+"_5_Guppy_4.0.11_prom.fastq.gz"

def get_sample(wildcards):
    return wildcards.sample

def get_out_folder(wildcards):
    return wildcards.sample+"/fastq_with_deletions"

rule get_pat_to_mat_assembly_insertions:
    input:
        bam="{sample}/pat_to_mat/{sample}.pat_subset_to_mat_subset.bam",
        bam_index="{sample}/pat_to_mat/{sample}.pat_subset_to_mat_subset.bam.bai",
    output:
        "{sample}/pat_to_mat/{sample}.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule run_convert_pat_to_mat_assembly_insertions_to_fasta:
    input:
        "{sample}/pat_to_mat/{sample}.insertions.tsv"
    output:
        temp("{sample}/pat_to_mat/{sample}.insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_assembly_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_pat_to_mat_assembly_insertion_annotation:
    input:
        tsv="{sample}/pat_to_mat/{sample}.insertions.tsv",
        tab="{sample}/pat_to_mat/{sample}.insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/pat_to_mat/{sample}.insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_assembly_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"



rule get_pat_and_mat_ref_insertions:
    input:
        mat_bam="{sample}/assembly_subset/{sample}.mat_subset_to_ref.bam",
        mat_bam_index="{sample}/assembly_subset/{sample}.mat_subset_to_ref.bam.bai",
        pat_bam="{sample}/assembly_subset/{sample}.pat_subset_to_ref.bam",
        pat_bam_index="{sample}/assembly_subset/{sample}.pat_subset_to_ref.bam.bai"
    output:
        mat="{sample}/assembly_subset/{sample}.mat.insertions.tsv",
        pat="{sample}/assembly_subset/{sample}.pat.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        """
        python {params.candidate_insertion_script} --bam {input.mat_bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output.mat}
        python {params.candidate_insertion_script} --bam {input.pat_bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output.pat}
        """

rule run_convert_pat_and_mat_ref_insertions_to_fasta:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.insertions.tsv",
        pat="{sample}/assembly_subset/{sample}.pat.insertions.tsv"
    output:
        mat=temp("{sample}/assembly_subset/{sample}.mat.insertions.fa"),
        pat=temp("{sample}/assembly_subset/{sample}.pat.insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_assembly_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        """
        python {params.candidate_insertion_conversion_script} --input {input.mat} > {output.mat}
        python {params.candidate_insertion_conversion_script} --input {input.pat} > {output.pat}
        """

rule run_pat_and_mat_ref_assembly_insertion_annotation:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.insertions.tsv",
        pat="{sample}/assembly_subset/{sample}.pat.insertions.tsv",
        mat_tab="{sample}/assembly_subset/{sample}.mat.insertions.mapped_to_repbase.last.tab",
        pat_tab="{sample}/assembly_subset/{sample}.pat.insertions.mapped_to_repbase.last.tab"
    output:
        mat="{sample}/assembly_subset/{sample}.mat.insertions.repbase_annotated.tsv",
        pat="{sample}/assembly_subset/{sample}.pat.insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_assembly_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        """
        python {params.candidate_insertion_annotation_script} --input {input.mat} --last {input.mat_tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output.mat}
        python {params.candidate_insertion_annotation_script} --input {input.pat} --last {input.pat_tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output.pat}
        """




rule get_mat_to_pat_assembly_insertions:
    input:
        bam="{sample}/mat_to_pat/{sample}.mat_subset_to_pat_subset.bam",
        bam_index="{sample}/mat_to_pat/{sample}.mat_subset_to_pat_subset.bam.bai",
    output:
        "{sample}/mat_to_pat/{sample}.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule run_convert_mat_to_pat_assembly_insertions_to_fasta:
    input:
        "{sample}/mat_to_pat/{sample}.insertions.tsv"
    output:
        temp("{sample}/mat_to_pat/{sample}.insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_assembly_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_mat_to_pat_assembly_insertion_annotation:
    input:
        tsv="{sample}/mat_to_pat/{sample}.insertions.tsv",
        tab="{sample}/mat_to_pat/{sample}.insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/mat_to_pat/{sample}.insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_assembly_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"



rule run_get_filtered_first_read_bams:
    input:
        mat_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.bam",
        mat_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.bam.bai",
        pat_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.bam",
        pat_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.bam.bai"
    output:
        mat_bam=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.1.bam"),
        mat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.merged_reads.1.txt"),
        pat_bam=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.1.bam"),
        pat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.merged_reads.1.txt")
    threads: 1
    params:
        merge_alignments_script = srcdir("../scripts/merge_split_alignments.py"),
        input_fastq=get_read_first_fastq,
        memory_per_thread="200G",
        tmp="{sample}/reads_assembly_subset_mapped/{sample}.tmp_output.1.bam"
    shell:
        """
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.mat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.mat_merged_reads}
        samtools sort {params.tmp} -T {output.mat_bam}.tmp -o {output.mat_bam}
        rm {params.tmp}
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.pat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.pat_merged_reads}
        samtools sort {params.tmp} -T {output.pat_bam}.tmp -o {output.pat_bam}
        rm {params.tmp}
        """

rule run_get_filtered_third_read_bams:
    input:
        mat_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.bam",
        mat_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.bam.bai",
        pat_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.bam",
        pat_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.bam.bai"
    output:
        mat_bam=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.3.bam"),
        mat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.merged_reads.3.txt"),
        pat_bam=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.3.bam"),
        pat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.merged_reads.3.txt")
    threads: 1
    params:
        merge_alignments_script = srcdir("../scripts/merge_split_alignments.py"),
        input_fastq=get_read_third_fastq,
        memory_per_thread="200G",
        tmp="{sample}/reads_assembly_subset_mapped/{sample}.tmp_output.3.bam"
    shell:
        """
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.mat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.mat_merged_reads}
        samtools sort {params.tmp} -T {output.mat_bam}.tmp -o {output.mat_bam}
        rm {params.tmp}
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.pat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.pat_merged_reads}
        samtools sort {params.tmp} -T {output.pat_bam}.tmp -o {output.pat_bam}
        rm {params.tmp}
        """

rule run_get_filtered_fifth_read_bams:
    input:
        mat_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.bam",
        mat_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.bam.bai",
        pat_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.bam",
        pat_bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.bam.bai"
    output:
        mat_bam=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.5.bam"),
        mat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.merged_reads.5.txt"),
        pat_bam=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.5.bam"),
        pat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.merged_reads.5.txt")
    threads: 1
    params:
        merge_alignments_script = srcdir("../scripts/merge_split_alignments.py"),
        input_fastq=get_read_fifth_fastq,
        memory_per_thread="200G",
        tmp="{sample}/reads_assembly_subset_mapped/{sample}.tmp_output.5.bam"
    shell:
        """
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.mat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.mat_merged_reads}
        samtools sort {params.tmp} -T {output.mat_bam}.tmp -o {output.mat_bam}
        rm {params.tmp}
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.pat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.pat_merged_reads}
        samtools sort {params.tmp} -T {output.pat_bam}.tmp -o {output.pat_bam}
        rm {params.tmp}
        """

rule get_mat_first_insertions:
    input:
        bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.1.bam",
        bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.1.bam.bai"
    output:
        "{sample}/reads_assembly_subset_mapped/{sample}.reads.mat.1.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule get_mat_third_insertions:
    input:
        bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.3.bam",
        bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.3.bam.bai"
    output:
        "{sample}/reads_assembly_subset_mapped/{sample}.reads.mat.3.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule get_mat_fifth_insertions:
    input:
        bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.5.bam",
        bam_index="{sample}/reads_assembly_subset_mapped/{sample}.mat.filtered.sorted.5.bam.bai",
    output:
        "{sample}/reads_assembly_subset_mapped/{sample}.reads.mat.5.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule get_pat_first_insertions:
    input:
        bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.1.bam",
        bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.1.bam.bai"
    output:
        "{sample}/reads_assembly_subset_mapped/{sample}.reads.pat.1.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule get_pat_third_insertions:
    input:
        bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.3.bam",
        bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.3.bam.bai"
    output:
        "{sample}/reads_assembly_subset_mapped/{sample}.reads.pat.3.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule get_pat_fifth_insertions:
    input:
        bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.5.bam",
        bam_index="{sample}/reads_assembly_subset_mapped/{sample}.pat.filtered.sorted.5.bam.bai",
    output:
        "{sample}/reads_assembly_subset_mapped/{sample}.reads.pat.5.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"


rule get_del_assemblies:
    input:
        mat_contigs="{sample}/assembly_subset/{sample}.mat.subset.fa",
        pat_contigs="{sample}/assembly_subset/{sample}.pat.subset.fa",
        mat_inserts="{sample}/assembly_subset/{sample}.mat.insertions.repbase_annotated.tsv",
        pat_inserts="{sample}/assembly_subset/{sample}.pat.insertions.repbase_annotated.tsv",
    output:
        mat_contigs="{sample}/assembly_subset/{sample}.mat.subset.deletions.fa",
        pat_contigs="{sample}/assembly_subset/{sample}.pat.subset.deletions.fa"
    params:
        script = srcdir("../scripts/delete_assembly_inserts.py"),
        memory_per_thread="16G"
    threads: 1
    shell:
        "python {params.script} --mat-ref-inserts {input.mat_inserts} --pat-ref-inserts {input.pat_inserts} --mat-contigs {input.mat_contigs} --pat-contigs {input.pat_contigs} --out-mat-contigs {output.mat_contigs} --out-pat-contigs {output.pat_contigs}"

rule zip_fastq:
    input:
        "{prefix}.fastq"
    output:
        "{prefix}.fastq.gz"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "bgzip -i {input}"

rule get_deletion_fastqs:
    input:
        mat_inserts="{sample}/assembly_subset/{sample}.mat.insertions.repbase_annotated.tsv",
        pat_inserts="{sample}/assembly_subset/{sample}.pat.insertions.repbase_annotated.tsv",
        reads_mat_1_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.bam",
        reads_mat_3_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.bam",
        reads_mat_5_bam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.bam",
        reads_pat_1_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.bam",
        reads_pat_3_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.bam",
        reads_pat_5_bam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.bam",
        reads_mat_1_bam_bai="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.bam.bai",
        reads_mat_3_bam_bai="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.bam.bai",
        reads_mat_5_bam_bai="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.bam.bai",
        reads_pat_1_bam_bai="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.bam.bai",
        reads_pat_3_bam_bai="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.bam.bai",
        reads_pat_5_bam_bai="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.bam.bai"
    output:
        tsv="{sample}/fastq_with_deletions/{sample}.annotation.tsv",
        f1="{sample}/fastq_with_deletions/Sample1.fastq",
        f3="{sample}/fastq_with_deletions/Sample3.fastq",
        f5="{sample}/fastq_with_deletions/Sample5.fastq",
        reads="{sample}/fastq_with_deletions/{sample}.reads_to_ignore.txt"
    threads: 1
    params:
        script=srcdir("../scripts/get_deletion_fastqs.py"),
        memory_per_thread="64G",
        first_fastq=get_read_first_fastq,
        third_fastq=get_read_third_fastq,
        fifth_fastq=get_read_fifth_fastq, 
        sample=get_sample,
        output_folder=get_out_folder
    shell:
        """
        python {params.script} --output-fasta-folder {params.output_folder} --reads-to-ignore {output.reads} --input-fastq-1 {params.first_fastq} --input-fastq-3 {params.third_fastq} --input-fastq-5 {params.fifth_fastq} --sample {params.sample} --mat-ref-inserts {input.mat_inserts} --pat-ref-inserts {input.pat_inserts} --mat-bam-1 {input.reads_mat_1_bam} --pat-bam-1 {input.reads_pat_1_bam} --mat-bam-3 {input.reads_mat_3_bam} --pat-bam-3 {input.reads_pat_3_bam} --mat-bam-5 {input.reads_mat_5_bam} --pat-bam-5 {input.reads_pat_5_bam} --centromeres {config[centromere_filter]} > {output.tsv}
        """
