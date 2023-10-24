import os
import glob

rule all_spike_in_bams:
    input:
        expand("HG{s}/HG{s}_{l}.spike_in.{r}.bam", s=config["samples"],l=config["read_lens"],r=config["reps"])

rule all_metrics:
    input:
        expand("simulation_results_HG{s}_{l}.txt" ,s=config["samples"],l=config["read_lens"])

include: "rules/simulation.smk"

