#
# Mapping related rules
#

def get_sample(wildcards):
    return wildcards.sample

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_base_dir(wildcards):
    return config["base_dir"]

def get_reference_base(wildcards):
    return config["reference_all"]

def get_reference_chr20(wildcards):
    return config["reference_chr20"]

def get_assembly_gfa(wildcards):
    return config["assembly_folder"]+"/"+wildcards.sample+"/assembly/"+wildcards.sample+".raw_unitig.gfa"

def get_first_fastq(wildcards):
    return wildcards.sample+"/fastq/"+wildcards.sample+"_1_Guppy_4.0.11_prom.fastq.gz"

def get_third_fastq(wildcards):
    return wildcards.sample+"/fastq/"+wildcards.sample+"_3_Guppy_4.0.11_prom.fastq.gz"

def get_control_fastq(wildcards):
    return wildcards.sample+"/fastq/control.fastq.gz"

def get_test_fastq(wildcards):
    return wildcards.sample+"/fastq/test.fastq.gz"

def get_fifth_fastq(wildcards):
    return wildcards.sample+"/fastq/"+wildcards.sample+"_5_Guppy_4.0.11_prom.fastq.gz"

def get_maternal(wildcards):
    return config["assembly_folder"]+"/"+wildcards.sample+"/assembly/"+wildcards.sample+".mat.fa"

def get_paternal(wildcards):
    return config["assembly_folder"]+"/"+wildcards.sample+"/assembly/"+wildcards.sample+".pat.fa"

def get_meryl_db(wildcards):
    return "merylDB_"+wildcards.sample

def get_fastq(wildcards):
    return "../validation/"+wildcards.sample+"/Sample1/fastq/Sample1.fastq.gz"

def get_mat_fa(wildcards):
    return wildcards.other+"/"+wildcards.other+"_chr20.mat.assembly.fa"

def get_pat_fa(wildcards):
    return wildcards.other+"/"+wildcards.other+"_chr20.pat.assembly.fa"

def get_prefix_other(wildcards):
    return "tmp_"+wildcards.sample+"_"+wildcards.other

def get_prefix(wildcards):
    return "tmp_"+wildcards.sample+"_combined"

rule minimap_align_assembly_mat:
    output:
        mat="{sample}/{sample}.reads_mapped_to.{other}.mat.bam"
    params:
        fastq=get_fastq,
        mat_fa=get_mat_fa,
        pat_fa=get_pat_fa,
        memory_per_thread="10G",
        prefix=get_prefix_other
    threads: 20
    shell:
        """
        minimap2 -ax map-ont -t {threads} {params.mat_fa} {params.fastq} | samtools sort > {output.mat}
        """

rule minimap_align_assembly_pat:
    output:
        pat="{sample}/{sample}.reads_mapped_to.{other}.pat.bam"
    params:
        fastq=get_fastq,
        mat_fa=get_mat_fa,
        pat_fa=get_pat_fa,
        memory_per_thread="10G",
        prefix=get_prefix_other
    threads: 20
    shell:
        """
        minimap2 -ax map-ont -t {threads} {params.pat_fa} {params.fastq} | samtools sort > {output.pat}
        """

rule minimap_align_grch38:
    output:
        grch="{sample}/{sample}.reads_mapped_to.grch38.bam"
    params:
        fastq=get_fastq,
        grch38_fa=get_reference_chr20,
        combined_fa="combined.chr20.assembly.fa",
        memory_per_thread="10G",
        prefix=get_prefix
    threads: 20
    shell:
        """
        minimap2 -ax map-ont -t {threads} {params.grch38_fa} {params.fastq} | samtools sort > {output.grch}
        """

rule minimap_align_combined:
    output:
        combined="{sample}/{sample}.reads_mapped_to.combined.bam"
    params:
        fastq=get_fastq,
        combined_fa="combined.chr20.assembly.fa",
        memory_per_thread="10G",
        prefix=get_prefix
    threads: 20
    shell:
        """
        minimap2 -ax map-ont -t {threads} --split-prefix {params.prefix} {params.combined_fa} {params.fastq} | samtools sort > {output.combined}
        """

# Winnowmap Alignmnet

rule winnow_index:
    input:
        get_reference_base
    output:
        config["reference_all"]+"_repetitive_k19.txt"
    params:
        memory_per_thread="48G"
    threads: 1
    shell:
        """
        {config[meryl_dir]} count k=19 output merylDB {input}
        {config[meryl_dir]} print greater-than distinct=0.9998 merylDB > {output}
        rm -r merylDB
        """


rule get_assembly_fa:
    output:
        "{sample}/assembly_mapped/{sample}.fa"
    threads: 1
    params:
        gfa= get_assembly_gfa,
        memory_per_thread="4G"
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' {params.gfa} > {output}
        """


rule winnow_align_assembly:
    input:
        fa="{sample}/assembly_mapped/{sample}.fa",
        winnow_index=config["reference_all"]+"_repetitive_k19.txt"
    output:
        protected("{sample}/assembly_mapped/{sample}.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use= get_reference_base
    threads: 20
    shell:
        """
        {config[winnow_dir]} -W {input.winnow_index} -ax asm20  {params.ref_to_use} {input.fa} | samtools sort -o {output}
        """


rule winnow_align_maternal:
    input:
        winnow_index=config["reference_all"]+"_repetitive_k19.txt"
    output:
        protected("{sample}/assembly_mapped/{sample}.mat.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use=get_reference_base,
        fasta=get_maternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -W {input.winnow_index} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule winnow_align_paternal:
    input:
        winnow_index=config["reference_all"]+"_repetitive_k19.txt"
    output:
        protected("{sample}/assembly_mapped/{sample}.pat.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use=get_reference_base,
        fasta=get_paternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -W {input.winnow_index} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """


rule winnow_index_assembly:
    input:
        "{sample}/assembly_mapped/{sample}.fa"
    output:
        "{sample}/assembly_mapped/{sample}_repetitive_k15.txt"
    params:
        memory_per_thread="48G",
        meryl_db= get_meryl_db
    threads: 1
    shell:
        """
        {config[meryl_dir]} count k=15 output {params.meryl_db} {input}
        {config[meryl_dir]} print greater-than distinct=0.9998 {params.meryl_db} > {output}
        rm -r {params.meryl_db}
        """

rule winnow_index_mat_assembly:
    output:
        "{sample}/assembly_mapped/{sample}_mat_repetitive_k15.txt"
    params:
        memory_per_thread="48G",
        meryl_db= get_meryl_db,
        fasta=get_maternal
    threads: 1
    shell:
        """
        {config[meryl_dir]} count k=15 output mat_{params.meryl_db} {params.fasta}
        {config[meryl_dir]} print greater-than distinct=0.9998 mat_{params.meryl_db} > {output}
        rm -r mat_{params.meryl_db}
        """

rule winnow_index_pat_assembly:
    output:
        "{sample}/assembly_mapped/{sample}_pat_repetitive_k15.txt"
    params:
        memory_per_thread="48G",
        meryl_db= get_meryl_db,
        fasta=get_paternal
    threads: 1
    shell:
        """
        {config[meryl_dir]} count k=15 output pat_{params.meryl_db} {params.fasta}
        {config[meryl_dir]} print greater-than distinct=0.9998 pat_{params.meryl_db} > {output}
        rm -r pat_{params.meryl_db}
        """

rule winnow_align_reads_to_assembly_first:
    input:
        assembly_fa="{sample}/assembly_mapped/{sample}.fa",
        winnow_index="{sample}/assembly_mapped/{sample}_repetitive_k15.txt",
        first_reads_fastq= get_first_fastq
    output:
        first=protected("{sample}/reads_assembly_mapped/{sample}.sorted.1.paf.gz")
    params:
        memory_per_thread="10G"
    threads: 20
    shell:
        """
        {config[winnow_dir]} -cx map-ont -W {input.winnow_index} -t {threads} {input.assembly_fa} {input.first_reads_fastq} | gzip > {output.first}
        """

rule winnow_align_reads_to_mat_assembly_first:
    input:
        winnow_index="{sample}/assembly_mapped/{sample}_mat_repetitive_k15.txt",
        first_reads_fastq= get_first_fastq
    output:
        first=protected("{sample}/reads_assembly_mapped/{sample}.mat.sorted.1.paf.gz")
    params:
        memory_per_thread="10G",
        fasta=get_maternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -cx map-ont -W {input.winnow_index} -t {threads} {params.fasta} {input.first_reads_fastq} | gzip > {output.first}
        """


rule winnow_align_reads_to_pat_assembly_first:
    input:
        winnow_index="{sample}/assembly_mapped/{sample}_pat_repetitive_k15.txt",
        first_reads_fastq= get_first_fastq
    output:
        first=protected("{sample}/reads_assembly_mapped/{sample}.pat.sorted.1.paf.gz")
    params:
        memory_per_thread="10G",
        fasta=get_paternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -cx map-ont -W {input.winnow_index} -t {threads} {params.fasta} {input.first_reads_fastq} | gzip > {output.first}
        """


rule winnow_align_reads_to_assembly_third:
    input:
        assembly_fa="{sample}/assembly_mapped/{sample}.fa",
        winnow_index="{sample}/assembly_mapped/{sample}_repetitive_k15.txt",
        third_reads_fastq= get_third_fastq
    output:
        third=protected("{sample}/reads_assembly_mapped/{sample}.sorted.3.paf.gz")
    params:
        memory_per_thread="10G"
    threads: 20
    shell:
        """
        {config[winnow_dir]} -W {input.winnow_index} -cx map-ont -t {threads} {input.assembly_fa} {input.third_reads_fastq} | gzip > {output.third}
        """


rule winnow_align_reads_to_mat_assembly_third:
    input:
        winnow_index="{sample}/assembly_mapped/{sample}_mat_repetitive_k15.txt",
        third_reads_fastq= get_third_fastq
    output:
        third=protected("{sample}/reads_assembly_mapped/{sample}.mat.sorted.3.paf.gz")
    params:
        memory_per_thread="10G",
        fasta=get_maternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -cx map-ont -W {input.winnow_index} -t {threads} {params.fasta} {input.third_reads_fastq} | gzip > {output.third}
        """


rule winnow_align_reads_to_pat_assembly_third:
    input:
        winnow_index="{sample}/assembly_mapped/{sample}_pat_repetitive_k15.txt",
        third_reads_fastq= get_third_fastq
    output:
        third=protected("{sample}/reads_assembly_mapped/{sample}.pat.sorted.3.paf.gz")
    params:
        memory_per_thread="10G",
        fasta=get_paternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -cx map-ont -W {input.winnow_index} -t {threads} {params.fasta} {input.third_reads_fastq} | gzip > {output.third}
        """


rule winnow_align_reads_to_assembly_fifth:
    input:
        assembly_fa="{sample}/assembly_mapped/{sample}.fa",
        winnow_index="{sample}/assembly_mapped/{sample}_repetitive_k15.txt",
        fifth_reads_fastq= get_fifth_fastq
    output:
        fifth=protected("{sample}/reads_assembly_mapped/{sample}.sorted.5.paf.gz")
    params:
        memory_per_thread="10G"
    threads: 20
    shell:
        """
        {config[winnow_dir]} -W {input.winnow_index} -cx map-ont -t {threads} {input.assembly_fa} {input.fifth_reads_fastq} | gzip > {output.fifth}
        """


rule winnow_align_reads_to_mat_assembly_fifth:
    input:
        winnow_index="{sample}/assembly_mapped/{sample}_mat_repetitive_k15.txt",
        fifth_reads_fastq= get_fifth_fastq
    output:
        fifth=protected("{sample}/reads_assembly_mapped/{sample}.mat.sorted.5.paf.gz")
    params:
        memory_per_thread="10G",
        fasta=get_maternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -cx map-ont -W {input.winnow_index} -t {threads} {params.fasta} {input.fifth_reads_fastq} | gzip > {output.fifth}
        """


rule winnow_align_reads_to_pat_assembly_fifth:
    input:
        winnow_index="{sample}/assembly_mapped/{sample}_pat_repetitive_k15.txt",
        fifth_reads_fastq= get_fifth_fastq
    output:
        fifth=protected("{sample}/reads_assembly_mapped/{sample}.pat.sorted.5.paf.gz")
    params:
        memory_per_thread="10G",
        fasta=get_paternal
    threads: 20
    shell:
        """
        {config[winnow_dir]} -cx map-ont -W {input.winnow_index} -t {threads} {params.fasta} {input.fifth_reads_fastq} | gzip > {output.fifth}
        """


rule paf_to_bam:
    input:
        fifth="{sample}/reads_assembly_mapped/{sample}.sorted.5.paf.gz",
        third="{sample}/reads_assembly_mapped/{sample}.sorted.3.paf.gz"
    output:
        fifth=protected("{sample}/reads_assembly_mapped/{sample}.sorted.5.bam"),
        third=protected("{sample}/reads_assembly_mapped/{sample}.sorted.3.bam")
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/paf_to_sam.py"),
        fifth_reads_fastq= get_fifth_fastq,
        third_reads_fastq= get_third_fastq,
        tmp_loc="{sample}/reads_assembly_mapped/{sample}.tmp",
        tmp_3="{sample}/reads_assembly_mapped/{sample}.3.bam",
        tmp_5="{sample}/reads_assembly_mapped/{sample}.5.bam"
    threads: 1
    shell:
        """
        python {params.script} --input-fastq {params.third_reads_fastq} --input-paf {input.third} --output-bam {params.tmp_3}
        python {params.script} --input-fastq {params.fifth_reads_fastq} --input-paf {input.fifth} --output-bam {params.tmp_5}
        samtools sort {params.tmp_3} -T {params.tmp_loc} -o {output.third}
        samtools sort {params.tmp_5} -T {params.tmp_loc} -o {output.fifth}
        rm {params.tmp_3}
        rm {params.tmp_5}
        """

rule paf_to_bam_mat_pat:
    input:
        fifth_mat="{sample}/reads_assembly_mapped/{sample}.mat.sorted.5.paf.gz",
        third_mat="{sample}/reads_assembly_mapped/{sample}.mat.sorted.3.paf.gz",
        fifth_pat="{sample}/reads_assembly_mapped/{sample}.pat.sorted.5.paf.gz",
        third_pat="{sample}/reads_assembly_mapped/{sample}.pat.sorted.3.paf.gz"
    output:
        fifth_mat=protected("{sample}/reads_assembly_mapped/{sample}.mat.sorted.5.bam"),
        third_mat=protected("{sample}/reads_assembly_mapped/{sample}.mat.sorted.3.bam"),
        fifth_pat=protected("{sample}/reads_assembly_mapped/{sample}.pat.sorted.5.bam"),
        third_pat=protected("{sample}/reads_assembly_mapped/{sample}.pat.sorted.3.bam")
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/paf_to_sam.py"),
        fifth_reads_fastq= get_fifth_fastq,
        third_reads_fastq= get_third_fastq,
        tmp_loc="{sample}/reads_assembly_mapped/{sample}.tmp",
        tmp_3="{sample}/reads_assembly_mapped/{sample}.3.bam",
        tmp_5="{sample}/reads_assembly_mapped/{sample}.5.bam"
    threads: 1
    shell:
        """
        python {params.script} --input-fastq {params.third_reads_fastq} --input-paf {input.third_mat} --output-bam {params.tmp_3}
        samtools sort {params.tmp_3} -T {params.tmp_loc} -o {output.third_mat}
        rm {params.tmp_3}
        python {params.script} --input-fastq {params.third_reads_fastq} --input-paf {input.third_pat} --output-bam {params.tmp_3}
        samtools sort {params.tmp_3} -T {params.tmp_loc} -o {output.third_pat}
        rm {params.tmp_3}
        python {params.script} --input-fastq {params.fifth_reads_fastq} --input-paf {input.fifth_pat} --output-bam {params.tmp_5}
        samtools sort {params.tmp_5} -T {params.tmp_loc} -o {output.fifth_pat}
        rm {params.tmp_5}
        python {params.script} --input-fastq {params.fifth_reads_fastq} --input-paf {input.fifth_mat} --output-bam {params.tmp_5}
        samtools sort {params.tmp_5} -T {params.tmp_loc} -o {output.fifth_mat}
        rm {params.tmp_5}
        """

rule paf_to_bam_sample1_mat_pat:
    input:
        first_mat="{sample}/reads_assembly_mapped/{sample}.mat.sorted.1.paf.gz",
        first_pat="{sample}/reads_assembly_mapped/{sample}.pat.sorted.1.paf.gz"
    output:
        first_mat=protected("{sample}/reads_assembly_mapped/{sample}.mat.sorted.1.bam"),
        first_pat=protected("{sample}/reads_assembly_mapped/{sample}.pat.sorted.1.bam")
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/paf_to_sam.py"),
        first_reads_fastq= get_first_fastq,
        tmp_loc="{sample}/reads_assembly_mapped/{sample}.tmp",
        tmp_1="{sample}/reads_assembly_mapped/{sample}.1.bam"
    threads: 1
    shell:
        """
        python {params.script} --input-fastq {params.first_reads_fastq} --input-paf {input.first_pat} --output-bam {params.tmp_1}
        samtools sort {params.tmp_1} -T {params.tmp_loc} -o {output.first_pat}
        rm {params.tmp_1}
        python {params.script} --input-fastq {params.first_reads_fastq} --input-paf {input.first_mat} --output-bam {params.tmp_1}
        samtools sort {params.tmp_1} -T {params.tmp_loc} -o {output.first_mat}
        rm {params.tmp_1}
        """


rule map_pat_to_mat:
    output:
        "{sample}/pat_to_mat/{sample}.pat_to_mat.bam"
    params:
        mat=get_maternal,
        pat=get_paternal,
        memory_per_thread="10G"
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {params.mat} {params.pat} | samtools sort -o {output}
        """

rule get_subset_assemblies:
    output:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa"
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_aligned_asssembly.py"),
        chr_list=get_chrom_list,
        sample=get_sample,
        input_path=get_base_dir
    threads: 1
    shell:
        """
        python {params.script} --chrom-list {params.chr_list} --sample {params.sample} --input-path {params.input_path}
        """

rule map_pat_subset_to_mat_subset:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa"
    output:
        "{sample}/pat_to_mat/{sample}.pat_subset_to_mat_subset.bam"
    params:
        memory_per_thread="10G"
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {input.mat} {input.pat} | samtools sort -o {output}
        """

rule map_mat_subset_to_pat_subset:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa"
    output:
        "{sample}/mat_to_pat/{sample}.mat_subset_to_pat_subset.bam"
    params:
        memory_per_thread="10G"
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {input.mat} {input.pat} | samtools sort -o {output}
        """


rule map_mat_subset_to_ref:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa"
    output:
        "{sample}/assembly_subset/{sample}.mat_subset_to_ref.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_reference_base
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {params.ref_to_use} {input.mat} | samtools sort -o {output}
        """

rule map_pat_subset_to_ref:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa"
    output:
        "{sample}/assembly_subset/{sample}.pat_subset_to_ref.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_reference_base
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {params.ref_to_use} {input.pat} | samtools sort -o {output}
        """

rule map_reads_to_pat_subset_control:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa",
        first_reads_fastq= get_control_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.control.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.control.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.control.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.pat} {input.first_reads_fastq} | samtools sort -o {output.first}
        """

rule map_reads_to_pat_subset_test:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa",
        first_reads_fastq= get_test_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.test.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.test.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.test.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.pat} {input.first_reads_fastq} | samtools sort -o {output.first}
        """


rule map_reads_to_mat_subset_control:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        first_reads_fastq= get_control_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.control.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.control.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.control.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.mat} {input.first_reads_fastq} | samtools sort -o {output.first}
        """

rule map_reads_to_mat_subset_test:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        first_reads_fastq= get_test_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.test.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.test.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.test.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.mat} {input.first_reads_fastq} | samtools sort -o {output.first}
        """



rule map_reads_to_pat_subset_1:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa",
        first_reads_fastq= get_first_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.pat} {input.first_reads_fastq} | samtools sort -o {output.first}
        """


rule map_reads_to_pat_subset_3:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa",
        third_reads_fastq= get_third_fastq
    output:
        third=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.pat} {input.third_reads_fastq} | samtools sort -o {output.third}
        """


rule map_reads_to_pat_subset_5:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.fa",
        fifth_reads_fastq= get_fifth_fastq
    output:
        fifth=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.pat} {input.fifth_reads_fastq} | samtools sort -o {output.fifth}
        """


rule map_reads_to_mat_subset_1:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        first_reads_fastq= get_first_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.mat} {input.first_reads_fastq} | samtools sort -o {output.first}
        """


rule map_reads_to_mat_subset_3:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        third_reads_fastq= get_third_fastq
    output:
        third=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.mat} {input.third_reads_fastq} | samtools sort -o {output.third}
        """


rule map_reads_to_mat_subset_5:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.fa",
        fifth_reads_fastq= get_fifth_fastq
    output:
        fifth=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.mat} {input.fifth_reads_fastq} | samtools sort -o {output.fifth}
        """


rule map_first_only:
    input:
        first=get_first_fastq
    output:
        first=protected("{sample}/reads_ref_mapped/{sample}.1.bam")
    params:
        memory_per_thread="10G",
        ref=get_reference_base
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax map-ont {params.ref} {input.first} | samtools sort > {output.first}
        """

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



rule make_lastdb:
    input:
        config["repbase"]
    output:
        config["repbase"] + ".lastdb.suf"
    threads: 1
    params:
        memory_per_thread="16G"
    shell:
        "lastdb {config[repbase]}.lastdb {input}"


rule map_insertion_sequences_last:
    input:
        fa="{base}.fa",
        db=config["repbase"] + ".lastdb.suf"
    output:
        "{base}.mapped_to_repbase.last.tab"
    threads: 8
    params:
        memory_per_thread="12G"
    shell:
        "lastal -P {threads} -f tab -r1 -a1 -b1 {config[repbase]}.lastdb {input.fa} > {output}"



rule map_reads_to_pat_del_1:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.deletions.fa",
        first_reads_fastq= get_first_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.deletions.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.deletions.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.1.deletions.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -cx map-ont {input.pat} {input.first_reads_fastq} | gzip > {params.tmp_paf}
        python {params.script} --input-fastq {input.first_reads_fastq} --input-paf {params.tmp_paf} --output-bam {params.tmp_sam}
        samtools sort {params.tmp_sam} -o {output.first}
        rm {params.tmp_paf}
        rm {params.tmp_sam}
        """


rule map_reads_to_pat_del_3:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.deletions.fa",
        third_reads_fastq= get_third_fastq
    output:
        third=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.deletions.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.deletions.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.3.deletions.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -cx map-ont {input.pat} {input.third_reads_fastq} | gzip > {params.tmp_paf}
        python {params.script} --input-fastq {input.third_reads_fastq} --input-paf {params.tmp_paf} --output-bam {params.tmp_sam}
        samtools sort {params.tmp_sam} -o {output.third}
        rm {params.tmp_paf}
        rm {params.tmp_sam}
        """


rule map_reads_to_pat_del_5:
    input:
        pat="{sample}/assembly_subset/{sample}.pat.subset.deletions.fa",
        fifth_reads_fastq= get_fifth_fastq
    output:
        fifth=protected("{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.deletions.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.deletions.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.pat.sorted.5.deletions.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -cx map-ont {input.pat} {input.fifth_reads_fastq} | gzip > {params.tmp_paf}
        python {params.script} --input-fastq {input.fifth_reads_fastq} --input-paf {params.tmp_paf} --output-bam {params.tmp_sam}
        samtools sort {params.tmp_sam} -o {output.fifth}
        rm {params.tmp_paf}
        rm {params.tmp_sam}
        """


rule map_reads_to_mat_del_1:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.deletions.fa",
        first_reads_fastq= get_first_fastq
    output:
        first=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.deletions.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.deletions.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.1.deletions.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -cx map-ont {input.mat} {input.first_reads_fastq} | gzip > {params.tmp_paf}
        python {params.script} --input-fastq {input.first_reads_fastq} --input-paf {params.tmp_paf} --output-bam {params.tmp_sam}
        samtools sort {params.tmp_sam} -o {output.first}
        rm {params.tmp_paf}
        rm {params.tmp_sam}
        """


rule map_reads_to_mat_del_3:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.deletions.fa",
        third_reads_fastq= get_third_fastq
    output:
        third=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.deletions.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.deletions.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.3.deletions.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -cx map-ont {input.mat} {input.third_reads_fastq} | gzip > {params.tmp_paf}
        python {params.script} --input-fastq {input.third_reads_fastq} --input-paf {params.tmp_paf} --output-bam {params.tmp_sam}
        samtools sort {params.tmp_sam} -o {output.third}
        rm {params.tmp_paf}
        rm {params.tmp_sam}
        """


rule map_reads_to_mat_del_5:
    input:
        mat="{sample}/assembly_subset/{sample}.mat.subset.deletions.fa",
        fifth_reads_fastq= get_fifth_fastq
    output:
        fifth=protected("{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.deletions.bam")
    params:
        memory_per_thread="10G",
        script=srcdir("../scripts/paf_to_sam.py"),
        tmp_paf="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.deletions.paf.gz",
        tmp_sam="{sample}/reads_assembly_subset_mapped/{sample}.mat.sorted.5.deletions.sam"
    threads: 20
    shell:
        """
        minimap2 -t {threads} -cx map-ont {input.mat} {input.fifth_reads_fastq} | gzip > {params.tmp_paf}
        python {params.script} --input-fastq {input.fifth_reads_fastq} --input-paf {params.tmp_paf} --output-bam {params.tmp_sam}
        samtools sort {params.tmp_sam} -o {output.fifth}
        rm {params.tmp_paf}
        rm {params.tmp_sam}
        """