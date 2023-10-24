
# Delete sequence to simulate somatic insertions

def get_sample(wildcards):
    return "wildcards.sample

def get_fastq(wildcards):
    return config["fastq"]

def get_ref(wildcards):
    return config["reference"]

def max_depth(wildcards):
    return 500

def get_base_dir(wildcards):
    return config["base_dir"]

def get_repbase(wildcards):
    return config["repbase"]

def get_centromeres(wildcards):
    return config["centromeres"]

def get_telomeres(wildcards):
    return config["telomeres"]

def get_maternal(wildcards):
    return config["mat_fa"]

def get_paternal(wildcards):
    return config["pat_fa"]

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_chrom_lengths(wildcards):
    return config["chrom_lengths"]

def get_pbsim_model(wildcards):
    return config["pbsim_model"]

def get_realign_tsv(wildcards):
    return config["run_dir"]+"/"+wildcards.sample+"/"+wildcards.read_len+"/RTD-Diplolid/test/results/test.inserts.tsv"

def get_sd(wildcards):
    return str(0.2*int(wildcards.read_len))

def get_pipeline_tsv(wildcards):
    return config["run_dir"]+"/"+wildcards.sample+"/"+wildcards.read_len+"/RTD-Diploid/test/winnow_realign_read_analysis/test.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv"

## Index a bam
rule make_bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "samtools index {input}"



rule align_maternal:
    output:
        "{sample}.mat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_maternal
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule align_paternal:
    output:
        "{sample}.pat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_paternal
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule get_insert_seqs:
    input:
        mat_bam="{sample}.mat.sorted.bam",
        mat_bai="{sample}.mat.sorted.bam.bai",
        pat_bam="{sample}.pat.sorted.bam",
        pat_bai="{sample}.pat.sorted.bam.bai"
    output:
        tsv="{sample}/{sample}.insert_added_sequences.tsv", 
        mat="{sample}/{sample}.mat.assembly_subset.fa",
        pat="{sample}/{sample}.pat.assembly_subset.fa"
    threads: 1
    params:
        memory_per_thread="24G",
        repbase_seqs=get_repbase,
        mat_fasta=get_maternal,
        pat_fasta=get_paternal,
        centromeres=get_centromeres,
        chrom_lengths=get_chrom_lengths,
        script=srcdir("../scripts/get_insert_seqs.py")
    shell:
        """
        python {params.script} --input {params.repbase_seqs} --mat-bam {input.mat_bam} --pat-bam {input.pat_bam} --mat-fasta {params.mat_fasta} --pat-fasta {params.pat_fasta} --centromeres {params.centromeres} --total 500 --chrom-lengths {params.chrom_lengths} --output-folder {wildcards.sample} --output-prefix {wildcards.sample}_inserted_sequence --mat-out-fa {output.mat} --pat-out-fa {output.pat} > {output.tsv}
        """


rule simulate_seqs:
    input:
        "{sample}/{sample}.insert_added_sequences.tsv"
    output:
        tsv="{sample}/{sample}_{read_len}.simulated_fastq_list.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model,
        script=srcdir("../scripts/generate_pbsim_run.py")
    shell:
        """
        python {params.script} --input {input} --output-folder {wildcards.sample} --output-prefix {wildcards.sample}_inserted_sequence --output-tsv {output.tsv} --pbsim-model {params.pbsim_model} --pbsim-path {config[pbsim]} --read-len {wildcards.read_len} > run_pbsim_{wildcards.sample}_{wildcards.read_len}.sh
        chmod +x run_pbsim_{wildcards.sample}_{wildcards.read_len}.sh
        ./run_pbsim_{wildcards.sample}_{wildcards.read_len}.sh
        """

rule find_spanning_reads:
    input:
        tsv=ancient("{sample}/{sample}_{read_len}.simulated_fastq_list.tsv"),
        inserts=ancient("{sample}/{sample}.insert_added_sequences.tsv")
    output:
        fastq="{sample}/{sample}_{read_len}.all_spanning_reads.fastq",
        tsv="{sample}/{sample}_{read_len}.spanning_reads_list.tsv"
    threads: 5
    params:
        memory_per_thread="15G",
        script=srcdir("../scripts/get_spanning_reads.py")
    shell:
        """
        python {params.script} --input {input.tsv} --inserts {input.inserts} --fastq {output.fastq} --tsv {output.tsv}
        """

rule simulate_fastq:
    input: 
        mat="{sample}/{sample}.mat.assembly_subset.fa",
        pat="{sample}/{sample}.pat.assembly_subset.fa"
    output:
        mat=temp("{sample}/{sample}_{read_len}.mat_full.{rep}.fastq"),
        pat=temp("{sample}/{sample}_{read_len}.pat_full.{rep}.fastq")
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model,
        sd=get_sd
    shell:
        """
        {config[pbsim]}/src/pbsim {input.mat} --prefix {wildcards.sample}_{wildcards.read_len}.mat_full.pbsim.{wildcards.rep} --difference-ratio 23:31:46 --seed {wildcards.rep} --length-min 1000 --depth 20 --hmm_model {config[pbsim]}{params.pbsim_model} --length-mean {wildcards.read_len} --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98 --length-sd {params.sd}
        rm {wildcards.sample}_{wildcards.read_len}.mat_full.pbsim.{wildcards.rep}*.maf
        rm {wildcards.sample}_{wildcards.read_len}.mat_full.pbsim.{wildcards.rep}*.ref
        cat {wildcards.sample}_{wildcards.read_len}.mat_full.pbsim.{wildcards.rep}* >> {wildcards.sample}/{wildcards.sample}_{wildcards.read_len}.mat_full.{wildcards.rep}.fastq
        rm {wildcards.sample}_{wildcards.read_len}.mat_full.pbsim.{wildcards.rep}*
        {config[pbsim]}/src/pbsim {input.pat} --prefix {wildcards.sample}_{wildcards.read_len}.pat_full.pbsim.{wildcards.rep} --difference-ratio 23:31:46 --seed {wildcards.rep} --length-min 1000 --depth 20 --hmm_model {config[pbsim]}{params.pbsim_model} --length-mean {wildcards.read_len} --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98 --length-sd {params.sd}
        rm {wildcards.sample}_{wildcards.read_len}.pat_full.pbsim.{wildcards.rep}*.maf
        rm {wildcards.sample}_{wildcards.read_len}.pat_full.pbsim.{wildcards.rep}*.ref
        cat {wildcards.sample}_{wildcards.read_len}.pat_full.pbsim.{wildcards.rep}* >> {wildcards.sample}/{wildcards.sample}_{wildcards.read_len}.pat_full.{wildcards.rep}.fastq
        rm {wildcards.sample}_{wildcards.read_len}.pat_full.pbsim.{wildcards.rep}*
        """

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


rule spike_fastqs:
    input:
        mat="{sample}/{sample}_{read_len}.mat_full.{rep}.fastq.gz",
        pat="{sample}/{sample}_{read_len}.pat_full.{rep}.fastq.gz",
        fastq="{sample}/{sample}_{read_len}.all_spanning_reads.fastq",
        tsv=ancient("{sample}/{sample}_{read_len}.spanning_reads_list.tsv")
    output:
        fastq="{sample}/{sample}_{read_len}.spike_in.{rep}.fastq",
        tsv="{sample}/{sample}_{read_len}.spike_in.{rep}.insert_list.txt"
    params:
        memory_per_thread="36G",
        script=srcdir("../scripts/spike_in_fastqs.py")
    shell:
        """
        python {params.script} --tsv {input.tsv} --inserts {input.fastq} --mat {input.mat} --pat {input.pat} --output-fastq {output.fastq} --rep {wildcards.rep} > {output.tsv}
        """


rule map_fastqs:
    input:
        fastq="{sample}/{sample}_{read_len}.spike_in.{rep}.fastq.gz"
    output:
        bam="{sample}/{sample}_{read_len}.spike_in.{rep}.bam"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {input.fastq} | samtools sort > {output.bam}
        """

rule summarize_data_spikein:
    output:
        tsv="simulation_results_{sample}_{read_len}.txt",
    threads: 1
    params:
        memory_per_thread="64G",
        somrit_test_script = srcdir("../scripts/get_simulation_metrics.py"),
        realign_tsv=get_realign_tsv,
        pipeline_tsv=get_pipeline_tsv,
        inserts_tsv="{sample}/{sample}.insert_added_sequences.tsv",
        spike_in_reads="{sample}/{sample}_{read_len}.spike_in.1.insert_list.txt"
    shell:
        """
        python {params.somrit_test_script} --inserts-tsv {params.inserts_tsv} --spiked-reads {params.spike_in_reads} --somrit {params.realign_tsv} --pipeline {params.pipeline_tsv} > {output.tsv}
        """

