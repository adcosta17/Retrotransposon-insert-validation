#
# Mapping related rules
#

def get_reference_base(wildcards):
    return config["reference_all"]

rule sniffles:
    input:
        bam = "{sample}/winnow_mapped/{sample}.sorted.bam",
        bai = "{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.sniffles.bedpe")
    params:
        memory_per_thread="10G",
        run_time="1:0:0:0",
        min_reads = "1",
        sniffles_dir = config["sniffles_dir"]
    threads: 10
    shell:
        "{params.sniffles_dir} -s {params.min_reads} --max_num_splits 1 -t {threads} -m {input.bam} -b {output} -n -1"