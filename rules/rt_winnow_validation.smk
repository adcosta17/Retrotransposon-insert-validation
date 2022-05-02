
# Validation rules

rule get_assembly_insertions:
    input:
        mat_bam="{sample}/assembly_mapped/{sample}.mat.sorted.bam",
        mat_bam_index="{sample}/assembly_mapped/{sample}.mat.sorted.bam.bai",
        pat_bam="{sample}/assembly_mapped/{sample}.pat.sorted.bam",
        pat_bam_index="{sample}/assembly_mapped/{sample}.pat.sorted.bam.bai"
    output:
        "{sample}/assembly_analysis/{sample}.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        "python {params.candidate_insertion_script} --mat-bam {input.mat_bam} --pat-bam {input.pat_bam} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"

rule run_convert_assembly_insertions_to_fasta_winnow:
    input:
        "{sample}/assembly_analysis/{sample}.insertions.tsv"
    output:
        temp("{sample}/assembly_analysis/{sample}.insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_assembly_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_assembly_insertion_annotation_winnow:
    input:
        tsv="{sample}/assembly_analysis/{sample}.insertions.tsv",
        tab="{sample}/assembly_analysis/{sample}.insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/assembly_analysis/{sample}.insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_assembly_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"


rule get_multi_sample_assembly:
    input:
        expand("{sample}/assembly_analysis/{sample}.insertions.repbase_annotated.tsv", sample=config["samples"]),
    output:
        "combined_assembly/combined_multi_sample.tsv"
    threads: 1
    params:
        gen_script = srcdir("../scripts/get_combined_assembly_inserts.py"),
        memory_per_thread="48G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} > {output}
        """


rule get_insert_supporting_reads:
    input:
        tsv="{sample}/assembly_analysis/{sample}.insertions.repbase_annotated.tsv",
        paf="{sample}/reads_assembly_mapped/{sample}.sorted.1.paf.gz"
    output:
        "{sample}/reads_assembly_mapped/{sample}.insert_supporting_reads.tsv"
    threads: 1
    params:
        script = srcdir("../scripts/get_reads_supporting_insertions.py"),
        memory_per_thread="24G"
    shell:
        "python {params.script} --input-tsv {input.tsv} --input-paf {input.paf} > {output}"


rule extract_read_insertions:
    input:
        tsv="{sample}/assembly_analysis/{sample}.insertions.repbase_annotated.tsv",
        paf="{sample}/reads_assembly_mapped/{sample}.sorted.1.paf.gz"
    output:
        "{sample}/reads_assembly_mapped/{sample}.extracted_insertions.tsv"
    threads: 1
    params:
        script = srcdir("../scripts/extract_insertions.py"),
        memory_per_thread="24G",
        fastq=get_first_fastq
    shell:
        "python {params.script} --input-tsv {input.tsv} --input-paf {input.paf} --input-fastq {params.fastq} > {output}"


rule extract_read_insertions_with_assembly_positions:
    input:
        tsv="{sample}/assembly_analysis/{sample}.insertions.repbase_annotated.tsv",
        mat_paf="{sample}/reads_assembly_mapped/{sample}.mat.sorted.1.paf.gz",
        pat_paf="{sample}/reads_assembly_mapped/{sample}.pat.sorted.1.paf.gz"
    output:
        "{sample}/reads_assembly_mapped/{sample}.extracted_insertions_with_positions.tsv"
    threads: 1
    params:
        script = srcdir("../scripts/extract_insertions_with_positions.py"),
        memory_per_thread="24G",
        fastq=get_first_fastq
    shell:
        "python {params.script} --input-tsv {input.tsv} --input-mat-paf {input.mat_paf} --input-pat-paf {input.pat_paf} --input-fastq {params.fastq} > {output}"


rule get_spliced_fastqs:
    input:
        supporting_reads = expand("{sample}/reads_assembly_mapped/{sample}.insert_supporting_reads.tsv", sample=config["samples"]),
        combined = "combined_assembly/combined_multi_sample.tsv"
    output:
        fastq_1 = "{sample}/fastq_inserted/{sample}_3.fastq.gz",
        fastq_2 = "{sample}/fastq_inserted/{sample}_5.fastq.gz",
        annotation = "{sample}/fastq_inserted/{sample}.insert_annotation.tsv"
    threads: 1
    params:
        script = srcdir("../scripts/generate_new_fastqs.py"),
        memory_per_thread="24G", 
        fastq_1 = "{sample}/fastq/{sample}_3_Guppy_4.0.11_prom.fastq.gz",
        fastq_2 = "{sample}/fastq/{sample}_5_Guppy_4.0.11_prom.fastq.gz",
        temp_fastq_1 = "{sample}/fastq_inserted/{sample}_3.fastq",
        temp_fastq_2 = "{sample}/fastq_inserted/{sample}_5.fastq",
    shell:
        """
        python {params.script} --input-fastq-1 {params.fastq_1} --input-fastq-2 {params.fastq_2} --output-fastq-1 {params.temp_fastq_1} --output-fastq-2 {params.temp_fastq_2} --sample {wildcards.sample} --all-samples {config[samples_csv]} --centromeres {config[centromere_filter]} --combined-tsv {input.combined} > {output.annotation}
        bgzip -c {params.temp_fastq_1} > {output.fastq_1}
        bgzip -c {params.temp_fastq_2} > {output.fastq_2}
        rm {params.temp_fastq_1}
        rm {params.temp_fastq_2}
        """


rule add_insertions_to_fastqs:
    input:
        supporting_reads = "{sample}/reads_assembly_mapped/{sample}.extracted_insertions.tsv"
    output:
        fastq_1 = "{sample}/fastq_insert_added/{sample}_3.fastq.gz",
        fastq_2 = "{sample}/fastq_insert_added/{sample}_5.fastq.gz",
        annotation = "{sample}/fastq_insert_added/{sample}.insert_annotation.tsv"
    threads: 1
    params:
        script = srcdir("../scripts/spike_in_fastqs.py"),
        memory_per_thread="24G", 
        fastq_1 = "{sample}/fastq/{sample}_3_Guppy_4.0.11_prom.fastq.gz",
        fastq_2 = "{sample}/fastq/{sample}_5_Guppy_4.0.11_prom.fastq.gz",
        temp_fastq_1 = "{sample}/fastq_insert_added/{sample}_3.fastq",
        temp_fastq_2 = "{sample}/fastq_insert_added/{sample}_5.fastq",
    shell:
        """
        python {params.script} --input-fastq-1 {params.fastq_1} --input-fastq-2 {params.fastq_2} --output-fastq-1 {params.temp_fastq_1} --output-fastq-2 {params.temp_fastq_2} --sample {wildcards.sample} --centromeres {config[centromere_filter]} --inserts-tsv {input.supporting_reads} > {output.annotation}
        bgzip -c {params.temp_fastq_1} > {output.fastq_1}
        bgzip -c {params.temp_fastq_2} > {output.fastq_2}
        rm {params.temp_fastq_1}
        rm {params.temp_fastq_2}
        """


rule add_insertions_to_fastq_positions:
    input:
        supporting_reads = "{sample}/reads_assembly_mapped/{sample}.extracted_insertions_with_positions.tsv",
        mat_bam_3="{sample}/reads_assembly_mapped/{sample}.mat.sorted.3.bam",
        mat_bam_5="{sample}/reads_assembly_mapped/{sample}.mat.sorted.5.bam",
        mat_bam_3_index="{sample}/reads_assembly_mapped/{sample}.mat.sorted.3.bam.bai",
        mat_bam_5_index="{sample}/reads_assembly_mapped/{sample}.mat.sorted.5.bam.bai",
        pat_bam_3="{sample}/reads_assembly_mapped/{sample}.pat.sorted.3.bam",
        pat_bam_5="{sample}/reads_assembly_mapped/{sample}.pat.sorted.5.bam",
        pat_bam_3_index="{sample}/reads_assembly_mapped/{sample}.pat.sorted.3.bam.bai",
        pat_bam_5_index="{sample}/reads_assembly_mapped/{sample}.pat.sorted.5.bam.bai",
        mat_asm_bam="{sample}/assembly_mapped/{sample}.mat.sorted.bam",
        mat_asm_bam_index="{sample}/assembly_mapped/{sample}.mat.sorted.bam.bai",
        pat_asm_bam="{sample}/assembly_mapped/{sample}.pat.sorted.bam",
        pat_asm_bam_index="{sample}/assembly_mapped/{sample}.pat.sorted.bam.bai"
    output:
        fastq_1 = "{sample}/fastq_insert_added_position/{sample}_3.fastq.gz",
        fastq_2 = "{sample}/fastq_insert_added_position/{sample}_5.fastq.gz",
        annotation = "{sample}/fastq_insert_added_position/{sample}.insert_annotation.tsv"
    threads: 1
    params:
        script = srcdir("../scripts/spike_in_fastqs_by_position.py"),
        memory_per_thread="24G", 
        fastq_1 = "{sample}/fastq/{sample}_3_Guppy_4.0.11_prom.fastq.gz",
        fastq_2 = "{sample}/fastq/{sample}_5_Guppy_4.0.11_prom.fastq.gz",
        temp_fastq_1 = "{sample}/fastq_insert_added_position/{sample}_3.fastq",
        temp_fastq_2 = "{sample}/fastq_insert_added_position/{sample}_5.fastq",
        mat_assembly_fa=get_maternal,
        pat_assembly_fa=get_paternal
    shell:
        """
        python {params.script} --input-fastq-1 {params.fastq_1} --input-fastq-2 {params.fastq_2} --output-fastq-1 {params.temp_fastq_1} --output-fastq-2 {params.temp_fastq_2} --sample {wildcards.sample} --centromeres {config[centromere_filter]} --mat-input-bam-1 {input.mat_bam_3} --mat-input-bam-2 {input.mat_bam_5} --pat-input-bam-1 {input.pat_bam_3} --pat-input-bam-2 {input.pat_bam_5} --pat-assembly-bam {input.pat_asm_bam} --mat-assembly-bam {input.mat_asm_bam} --inserts-tsv {input.supporting_reads} --mat-assembly-fa {params.mat_assembly_fa} --pat-assembly-fa {params.pat_assembly_fa} > {output.annotation}
        bgzip -c {params.temp_fastq_1} > {output.fastq_1}
        bgzip -c {params.temp_fastq_2} > {output.fastq_2}
        rm {params.temp_fastq_1}
        rm {params.temp_fastq_2}
        """
