# Retrotransposon-insert-validation

A set of snakemake rules and scripts for validating the calls generated by the [Retrotransposon Insertion Detection Pipeline](https://github.com/adcosta17/Retrotransposon-insert-detection)

## Dependencies
- pysam
- bgzip
- snakemake
- minimap2
- winnowmap2
- lastal
- samtools
- longshot
- whatshap

## Overview

This repo contains a set of snakemake rules and scripts that are used to generate the data used in, and evaluate the results of experiments designed to assess the performance of the [Retrotransposon Insertion Detection Pipeline](https://github.com/adcosta17/Retrotransposon-insert-detection).

### Deletion Based Validation

A deletion based validation approach is used to assess the precision and recall of the pipeline using PacBio Hifi based assemblies and Nanopore long reads.
We downloaded the PacBio Hifi based diploid assemblies and Nanopore reads from the [HPRC](https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0) for samples: [HG00438](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/HG00438/), [HG00621](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/HG00621/), [HG00673](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/HG00673/), [HG00735](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/HG00735/), and [HG00741](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/HG00741/).
Please see the methods section of the manuscript for a full breakdown of the validation methods. 

#### Prerequisites & Setup

To generate modified fastqs the snakemake pipeline requires a folder for each sample. Within each sample folder a subfolder called `fastq` that contains the original fastqs is also required. The fastqs should be compressed using bgzip for pysam to be able to read them. Additionally a second subfolder called `assembly` should be generated that contains the maternal and paternal contigs downloaded from the HPRC. The maternal should be named `<sample>.mat.fa` and the paternal `<sample>.pat.fa`. 

Each HPRC sample listed contains three fastqs each. These can be downloaded directly to the sample's fastq folder. We can then generate deletion fastqs by listing each sample in the config and running the snakemake listed below:

```sh
# Generate Deletion fastqs for HPRC Samples
snakemake -s validation.smk --configfile project_config.yaml all_deletion_fastq
```

#### Results


Once the deletion fastqs have been generated, the main retrotransposon insertion detection pipeline can be run on the modified fastqs. The results of this pipeline can be compared to the annotation file generated when the deletion fastqs are created.

```sh
# Get TP, FP, and FN values for Other HPRC Samples deletion based validation runs
python scripts/get_combined_metrics_deletions.py --input-prefix <path/to/rt_insert_pipeline_run_folder> --id-list <samples,csv,list> --folder winnow_realign_read_analysis --bam-folder winnow_realign
```

### False Positive Rate Validation

In addition to the deletion validation we ran the pipeline on the original unmodified fastqs for each of the HPRC samples listed above as well as data from [HG0002](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC_PLUS/HG002/). For the HPRC samples listed above we used the same setup to run the pipeline, one control and two tests,  but did not make any modifications to the fastqs. For HG0002 we combined two downloaded fastqs to generate a coverage of ~ 100x and then randomly split this into smaller fastqs of roughly equal coverage, shown below, before running the pipeline on the fastq splits selecting one as a control and the other as tests.  

As all the samples run were assumed to be normal samples free of any somatic variation, any insertions detected in the test fastqs are false positives. 

#### HG0002

The fastqs downloaded from the HPRC are quite large and need to be downsampled to generate multiple fastqs of equal lower coverage.

```sh
# Generate 3 fastqs of 20x each from a larger fastq of 60x
python scripts/get_HG002_fastqs.py --input-fastq <input_HG0002.fastq> --output-prefix <output_folder> --coverage 20 --test-number 3 --total-coverage 60
```

Once downsampled fastqs have been generated we can create a folder for the sample at the requested coverage and  copy the downsampled fastqs to the fastq subfolder. 






