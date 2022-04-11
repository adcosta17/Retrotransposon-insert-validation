#
# Mapping related rules for short-reads
#

def get_reference_base(wildcards):
    return config["reference_all"]

def get_cram(wildcards):
    return wildcards.sample+"/illumina/"+wildcards.sample+".final.cram"

rule get_short_read_bams:
    input:
        get_reference_base
    output:
        temp("{sample}/illumina/{sample}.bam")
    threads: 1
    params:
        ref= get_reference_base,
        cram= get_cram,
        memory_per_thread="12G"
    shell:
        """
        samtools view -b {params.cram} -T {params.ref} > {output}
        """

rule get_downsampled_fastqs:
    input:
        p1="{sample}/illumina/{sample}.p1.fastq.gz",
        p2="{sample}/illumina/{sample}.p2.fastq.gz"
    output:
        p1="{sample}/illumina/{sample}.p1.0.5.fastq.gz",
        p2="{sample}/illumina/{sample}.p2.0.5.fastq.gz"
    threads: 1
    params:
        memory_per_thread="12G"
    shell:
        """
        seqtk sample -s100 {input.p1} 0.5 | gzip > {output.p1}
        seqtk sample -s100 {input.p2} 0.5 | gzip > {output.p2}
        """

rule abyss_assemble:
    output:
        "{sample}/illumina_assembly/assembly-scaffolds.fa"
    threads: 10
    params:
        memory_per_thread="40G",
        p1="{sample}/illumina/p1.fastq.gz",
        p2="{sample}/illumina/p2.fastq.gz"
    shell:
        """
        abyss-pe np={threads} name={wildcards.sample}/illumina/assembly k=96 in=\'{params.p1} {params.p2}\' B=50G H=4 kc=4 v=-v
        """
