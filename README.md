README
================
2023-05-04

# Quantifying Isoform Diversity with lrRNA-Seq in Mouse Brain

## Purpose

The purpose of this project is to analyze Oxford Nanopore RNA sequencing
from four wildtype mouse brain regions balanced for sex, to assay
isoform usage differences across brain region and sex.

## Dependencies

“Dependencies” Include a code chunk (if relevant) on what your project
depends on.  
Currently: nextflow, gffread

## Scripts

“Scripts” Include a “script tree” on how to reproduce your results with
descriptions for all scripts

    ## src
    ## ├── README
    ## ├── dtu_analyses
    ## │   ├── calculate_cpm.R
    ## │   ├── dtu_isoform_switching.Rmd
    ## │   ├── dtu_neuro_diseases.Rmd
    ## │   ├── dtu_region_others.Rmd
    ## │   ├── dtu_region_region.Rmd
    ## │   ├── dtu_region_sex.Rmd
    ## │   ├── functions.R
    ## │   ├── pca_eda.Rmd
    ## │   └── size_power.R
    ## └── preprocessing
    ##     ├── 20010631.c0158.err.txt
    ##     ├── 20010631.c0158.out.txt
    ##     ├── README
    ##     ├── create_isoform_fa.sh
    ##     ├── get_genome_annotations.sh
    ##     └── run_nanoseq.sh

### Reproducibility

HPC user session info: I used the short partition with 1-2 nodes and
80GB per node, which is generous. None of my scripts take longer than 12
hours so I would stay on the short partition, unless you want to have
one long session on another partition which is fine.

*Script 00: preprocessing/run_nanoseq.sh* - This script can be skipped
if you are working from counts files. Does not require a docker. You
need raw input fasta files, a genome file, and genome annotation file. I
used GENCODE mouse release M31.

**Start here**

*Script 01: dtu_analyses/calculate_cpm.R* - This script is very short
and depends on having transcript counts available either by downloading
or running the nanoseq pipeline. Can run in either docker. Takes less
than one minute.

*Script 02: preprocessing/get_genome_annotations.sh* - you need these
for gffread to run, you also need genome annotations to run the nanoseq
pipeline so frankly you should already have them somewhere. You do need
to get genome annotations if you didn’t run the nanoseq pipeline and are
just using counts, so I am keeping it as script 2. It took me 20 minutes
to download.

*Script 03: preprocessing/create_isoform_fa.sh* - this script depends on
script 02. It is a single command but needs to be run in a docker
(either RStudio docker I made has gffread) or on a local machine with
gffread.

*Script 04: dtu_analyses/pca_eda.* This script is for exploratory data
analysis and PCA. It depends on script 01. If data looks bad, do not
proceed to script 05, but it is technically independent. This script
takes 3 minutes to run.

*Script 05: dtu_analyses/dtu_region_region.R* - This script depends on
the outputs from script 01 which are read as an RDS file. This script
also depends on script 03 - the gffread script. Run in github.1.2 docker
so you have the most up-to-date version of the package. This script
takes about an hour to run.

## Authors

Emma Jones, TC Howton, Victoria Flanary, and Brittany Lasseigne.

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

## Funding

NHGRI R00HG009678 (PI: Lasseigne)  
Pittman Scholar (PI: Lasseigne)  
UAB Lab Startup funds (PI: Lasseigne)

## Acknowledgements

We acknowledge the members of the Lasseigne Lab for their thoughtful
feedback.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.
