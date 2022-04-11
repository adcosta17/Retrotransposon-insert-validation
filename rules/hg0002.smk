
# Validation rules

def get_read_fastq(wildcards):
    return wildcards.sample+"/fastq/"+wildcards.sample+"_1_Guppy_4.0.11_prom.fastq.gz"

def get_sample(wildcards):
    return wildcards.sample

def get_out_folder(wildcards):
    return wildcards.sample+"/fastq_with_deletions"

def get_fastq_folder(wildcards):
    return wildcards.sample+"/fastq/"

def get_fastq(wildcards):
    return wildcards.sample+"/fastq/HG0002."+wildcards.i+".fastq.gz"

def get_base_dir(wildcards):
    return config['base_dir']

def get_mat_fa(wildcards):
    return config['base_dir']+"/HG0002/assembly/HG0002.mat.fa"

def get_pat_fa(wildcards):
    return config['base_dir']+"/HG0002/assembly/HG0002.pat.fa"

def get_ref(wildcards):
    return config["reference_all"]

rule minimap_align_assembly_mat_hg00002:
    output:
        mat="HG0002/assembly_mapped/HG0002.mat.sorted.bam"
    params:
        ref=get_ref,
        mat_fa=get_mat_fa,
        memory_per_thread="10G"
    threads: 20
    shell:
        """
        minimap2 -ax asm20 -t {threads} {params.ref} {params.mat_fa} | samtools sort > {output.mat}
        """

rule minimap_align_assembly_pat_hg00002:
    output:
        pat="HG0002/assembly_mapped/HG0002.pat.sorted.bam"
    params:
        ref=get_ref,
        pat_fa=get_pat_fa,
        memory_per_thread="10G"
    threads: 20
    shell:
        """
        minimap2 -ax asm20 -t {threads} {params.ref} {params.pat_fa} | samtools sort > {output.pat}
        """

rule get_subset_assemblies_hg00002:
    input:
        mat="HG0002/assembly_mapped/HG0002.mat.sorted.bam",
        pat="HG0002/assembly_mapped/HG0002.pat.sorted.bam",
        mat_idx="HG0002/assembly_mapped/HG0002.mat.sorted.bam.bai",
        pat_idx="HG0002/assembly_mapped/HG0002.pat.sorted.bam.bai"
    output:
        mat="HG0002/assembly_subset/HG0002.mat.subset.fa",
        pat="HG0002/assembly_subset/HG0002.pat.subset.fa"
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_aligned_asssembly.py"),
        chr_list=get_chrom_list,
        input_path=get_base_dir
    threads: 1
    shell:
        """
        python {params.script} --chrom-list {params.chr_list} --sample HG0002 --input-path {params.input_path}
        """

rule map_mat_subset_to_ref_hg00002:
    input:
        mat="HG0002/assembly_subset/HG0002.mat.subset.fa"
    output:
        "HG0002/assembly_subset/HG0002.mat_subset_to_ref.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {params.ref_to_use} {input.mat} | samtools sort -o {output}
        """

rule map_pat_subset_to_ref_hg00002:
    input:
        pat="HG0002/assembly_subset/HG0002.pat.subset.fa"
    output:
        "HG0002/assembly_subset/HG0002.pat_subset_to_ref.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {params.ref_to_use} {input.pat} | samtools sort -o {output}
        """

rule map_reads_to_mat_hg00002:
    input:
        mat="HG0002/assembly_subset/HG0002.mat.subset.fa",
    output:
        bam=protected("{sample}/reads_assembly_subset_mapped/HG0002.mat.sorted.{i}.bam")
    params:
        memory_per_thread="10G",
        fastq=get_fastq
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.mat} {params.fastq} | samtools sort -o {output.bam}
        """

rule map_reads_to_pat_hg00002:
    input:
        pat="HG0002/assembly_subset/HG0002.pat.subset.fa",
    output:
        bam=protected("{sample}/reads_assembly_subset_mapped/HG0002.pat.sorted.{i}.bam")
    params:
        memory_per_thread="10G",
        fastq=get_fastq
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.pat} {params.fastq} | samtools sort -o {output.bam}
        """

rule get_pat_and_mat_ref_insertions_hg00002:
    input:
        mat_bam="HG0002/assembly_subset/HG0002.mat_subset_to_ref.bam",
        mat_bam_index="HG0002/assembly_subset/HG0002.mat_subset_to_ref.bam.bai",
        pat_bam="HG0002/assembly_subset/HG0002.pat_subset_to_ref.bam",
        pat_bam_index="HG0002/assembly_subset/HG0002.pat_subset_to_ref.bam.bai"
    output:
        mat="HG0002/assembly_subset/HG0002.mat.insertions.tsv",
        pat="HG0002/assembly_subset/HG0002.pat.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        """
        python {params.candidate_insertion_script} --bam {input.mat_bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output.mat}
        python {params.candidate_insertion_script} --bam {input.pat_bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output.pat}
        """

rule run_convert_pat_and_mat_ref_insertions_to_fasta_hg00002:
    input:
        mat="HG0002/assembly_subset/HG0002.mat.insertions.tsv",
        pat="HG0002/assembly_subset/HG0002.pat.insertions.tsv"
    output:
        mat=temp("HG0002/assembly_subset/HG0002.mat.insertions.fa"),
        pat=temp("HG0002/assembly_subset/HG0002.pat.insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_assembly_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        """
        python {params.candidate_insertion_conversion_script} --input {input.mat} > {output.mat}
        python {params.candidate_insertion_conversion_script} --input {input.pat} > {output.pat}
        """

rule run_pat_and_mat_ref_assembly_insertion_annotation_hg00002:
    input:
        mat="HG0002/assembly_subset/HG0002.mat.insertions.tsv",
        pat="HG0002/assembly_subset/HG0002.pat.insertions.tsv",
        mat_tab="HG0002/assembly_subset/HG0002.mat.insertions.mapped_to_repbase.last.tab",
        pat_tab="HG0002/assembly_subset/HG0002.pat.insertions.mapped_to_repbase.last.tab"
    output:
        mat="HG0002/assembly_subset/HG0002.mat.insertions.repbase_annotated.tsv",
        pat="HG0002/assembly_subset/HG0002.pat.insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_assembly_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        """
        python {params.candidate_insertion_annotation_script} --input {input.mat} --last {input.mat_tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output.mat}
        python {params.candidate_insertion_annotation_script} --input {input.pat} --last {input.pat_tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output.pat}
        """


rule run_get_filtered_read_bams_hg00002:
    input:
        mat_bam="{sample}/reads_assembly_subset_mapped/HG0002.mat.sorted.{i}.bam",
        mat_bam_index="{sample}/reads_assembly_subset_mapped/HG0002.mat.sorted.{i}.bam.bai",
        pat_bam="{sample}/reads_assembly_subset_mapped/HG0002.pat.sorted.{i}.bam",
        pat_bam_index="{sample}/reads_assembly_subset_mapped/HG0002.pat.sorted.{i}.bam.bai"
    output:
        mat_bam=protected("{sample}/reads_assembly_subset_mapped/HG0002.mat.filtered.sorted.{i}.bam"),
        mat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/HG0002}.mat.merged_reads.{i}.txt"),
        pat_bam=protected("{sample}/reads_assembly_subset_mapped/HG0002.pat.filtered.sorted.{i}.bam"),
        pat_merged_reads=protected("{sample}/reads_assembly_subset_mapped/HG0002.pat.merged_reads.{i}.txt")
    threads: 1
    params:
        merge_alignments_script = srcdir("../scripts/merge_split_alignments.py"),
        input_fastq=get_fastq,
        memory_per_thread="200G",
        tmp="{sample}/reads_assembly_subset_mapped/HG0002.tmp_output.{i}.bam"
    shell:
        """
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.mat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.mat_merged_reads}
        samtools sort {params.tmp} -T {output.mat_bam}.tmp -o {output.mat_bam}
        rm {params.tmp}
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.pat_bam} --indel-size {config[min_insertion_length]} --fastq {params.input_fastq} --bam-output {params.tmp} --merged {output.pat_merged_reads}
        samtools sort {params.tmp} -T {output.pat_bam}.tmp -o {output.pat_bam}
        rm {params.tmp}
        """

rule get_deletion_fastqs_hg00002:
    input:
        mat_inserts="HG0002/assembly_subset/HG0002.mat.insertions.repbase_annotated.tsv",
        pat_inserts="HG0002/assembly_subset/HG0002.pat.insertions.repbase_annotated.tsv",
        reads_mapped_mat=expand("{{sample}}/reads_assembly_subset_mapped/HG0002.mat.sorted.{i}.bam",i=config["counts"]),
        reads_mapped_pat=expand("{{sample}}/reads_assembly_subset_mapped/HG0002.pat.sorted.{i}.bam",i=config["counts"]),
        reads_mapped_mat_bai=expand("{{sample}}/reads_assembly_subset_mapped/HG0002.mat.sorted.{i}.bam.bai",i=config["counts"]),
        reads_mapped_pat_bai=expand("{{sample}}/reads_assembly_subset_mapped/HG0002.pat.sorted.{i}.bam.bai",i=config["counts"])
    output:
        tsv="{sample}/fastq_with_deletions/HG0002.annotation.tsv",
        fastqs=expand("{{sample}}/fastq_with_deletions/HG0002.{i}.fastq", i=config["counts"]),
        reads="{sample}/fastq_with_deletions/HG0002.reads_to_ignore.txt"
    threads: 1
    params:
        script=srcdir("../scripts/get_HG002_deletion_fastqs.py"),
        memory_per_thread="64G",
        fastq_folder=get_fastq_folder,
        output_folder=get_out_folder
    shell:
        """
        python {params.script} --output-fasta-folder {params.output_folder} --mat-ref-inserts {input.mat_inserts} --pat-ref-inserts {input.pat_inserts} --centromeres {config[centromere_filter]} --sample {wildcards.sample} --counts {config[total_counts]} --input-folder {params.fastq_folder} > {output.tsv}
        """
