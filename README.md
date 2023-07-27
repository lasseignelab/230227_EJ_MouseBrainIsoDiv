README
================
2023-07-27

# Quantifying Isoform-Level Diversity with lrRNA-Seq in WT Mouse Brain

## Authors

Emma Jones, TC Howton, Victoria Flanary, and Brittany Lasseigne.

## Purpose

The purpose of this project is to analyze Oxford Nanopore RNA sequencing
from four wildtype mouse brain regions balanced for sex, to assay
isoform usage differences across brain region and sex.

## Scripts

    ## src
    ## ├── README
    ## ├── de_analysis
    ## │   ├── 10_DESeq2_region_region.Rmd
    ## │   ├── 11_DESeq2_region_others.Rmd
    ## │   ├── 12_DESeq2_region_sex.Rmd
    ## │   ├── 13_incorporate_de_results.Rmd
    ## │   ├── 15_compare_DTU_DGE.Rmd
    ## │   └── de_functions.R
    ## ├── dtu_analyses
    ## │   ├── 01_calculate_cpm.R
    ## │   ├── 04_pca_eda.Rmd
    ## │   ├── 05_dtu_region_region.Rmd
    ## │   ├── 06_dtu_region_others.Rmd
    ## │   ├── 07_dtu_region_sex.Rmd
    ## │   ├── 08_dtu_neuro_diseases.Rmd
    ## │   ├── 09_dtu_isoform_switching.Rmd
    ## │   ├── 14_protein_domain_info.Rmd
    ## │   ├── compare_satuRn_DEXSeq_results.Rmd
    ## │   ├── functions.R
    ## │   ├── isoformswitchanalyzer_overflow.Rmd
    ## │   └── size_power.R
    ## ├── preprocessing
    ## │   ├── 00_run_nanoseq.sh
    ## │   ├── 02_get_genome_annotations.sh
    ## │   ├── 03_create_isoform_fa.sh
    ## │   ├── 16_dataset_overview_figure.Rmd
    ## │   ├── 20010631.c0158.err.txt
    ## │   ├── 20010631.c0158.out.txt
    ## │   └── README
    ## └── worm_analysis
    ##     ├── 17_dtu_whole_worm.Rmd
    ##     ├── 18_compare_worm_mouse.Rmd
    ##     └── salmon_quant.sh

### Reproducibility

HPC user session info: I used the short partition with 1-2 nodes and
80GB per node, which is generous. None of my scripts take longer than 12
hours so I would stay on the short partition, unless you want to have
one long session on another partition which is fine.

***Script 00: preprocessing/00_run_nanoseq.sh*** - This script can be
skipped if you are working from counts files. Does not require a docker.
You need raw input fasta files, a genome file, and genome annotation
file. I used GENCODE mouse release M31.

#### Start here

***Script 01: dtu_analyses/01_calculate_cpm.R*** - This script is very
short and depends on having transcript counts available either by
downloading or running the nanoseq pipeline. Takes less than one minute.

***Script 02: preprocessing/02_get_genome_annotations.sh*** - You need
genome annotations available for gffread to run, but you also need
genome annotations to run the nanoseq pipeline so frankly you should
already have them somewhere and could just move them. You do need to get
genome annotations if you didn’t run the nanoseq pipeline and are just
using counts, so I am keeping it as script 2. It took me 20 minutes to
download.

***Script 03: preprocessing/03_create_isoform_fa.sh*** - This script
depends on script 02. It is a single command but needs to be run in a
docker (either RStudio docker I made has gffread) or on a local machine
with gffread.

***Script 04: dtu_analyses/04_pca_eda.Rmd*** - This script is for
exploratory data analysis and PCA. It depends on script 01. If data
looks bad, do not proceed to script 05, but it is technically
independent. This script takes 3 minutes to run.

***Script 05: dtu_analyses/05_dtu_region_region.Rmd*** - This script
depends on the outputs from script 01 which are read as an RDS file.
This script also depends on script 03 - the gffread script and the
functions script. Run in github.1.2 docker so you have the most
up-to-date version of the package. This script takes about 10 minutes to
run.

***Script 06: dtu_analyses/06_dtu_region_others.Rmd*** - This script is
meant for DTU analysis comparing one brain region to all others. This
script depends on scripts 01-03 and the functions script. Run this in
the github.1.2 docker. Takes about 10 minutes to run.

***Script 07: dtu_analyses/07_dtu_region_sex.Rmd*** - This script is
meant for DTU analysis across sexes within brain regions. This script
depends on scripts 01-03 and functions script. Run this in the
github.1.2 docker. Takes about 20 minutes to run (I think).

***Script 08: dtu_analyses/08_dtu_neuro_diseases.Rmd*** - This script is
for comparing brain-region-specific DTU genes to known Alzheimer’s
disease, psychiatric disorder, and CPAM case genes. It is dependent on
scripts 01 - 06 and the functions script. It takes 4 minutes to run.

***Script 09: dtu_analyses/09_dtu_isoform_switching.Rmd*** - This script
is for adding open reading frames to the existing switchlist objects. It
is also for plotting and saving switch plots, which show the significant
isoform switching events. It is dependent on scripts 01-07 and the
functions script. This script also takes 36 minutes to run.

***Script 10: de_analysis/10_DESeq2_region_region.Rmd*** - This script
is dependent on script 00/01 or having gene and transcript level count
data with metadata in your /data/ directory and the de functions script.
The purpose of this script is to run DESeq2 across brain regions at the
gene and transcript level. It takes less than 5 minutes to run.

***Script 11: de_analysis/11_DESeq2_region_others.Rmd*** - This script
is dependent on script 00/01 or having gene and transcript level count
data with metadata in your /data/ directory and the de functions script.
The purpose of this script is to run DESeq2 for each brain region
compared to an aggregate of other brain regions at the gene and
transcript level. It takes about 10 minutes to run.

***Script 12: de_analysis/12_DESeq2_region_sex.Rmd*** - This script is
dependent on script 00/01 or having gene and transcript level count data
with metadata in your /data/ directory and the de functions script. The
purpose of this script is to compare expression at the gene and
transcript level across sexes within brain regions. It takes less than
10 minutes to run.

***Script 13: de_analysis/13_incorporate_de_results.Rmd*** - This script
is fully dependent on scripts 01-09 and the de_analysis DESeq2 scripts
10-12 and the de functions script. The purpose of this script is to wrap
in all DESeq2 significance values into my isoformSwitchAnalyzeR objects
and subsequently, plots. This enables us to plot them all with
significance values for DGE, DTE, and DTU. This script takes 12 minutes
to run.

***Script 14: dtu_analyses/14_protein_domain_info.Rmd*** - This script
depends on scripts 00-13, and also includes a shell script for running
perl code. The purpose of this script is to extract nucleotide and amino
acid sequences and run pfam to annotate the protein domains. NEEDS TO
RUN IN PFAM DOCKER!

***Script 15: de_analysis/15_compare_DTU_DGE.Rmd*** - This script
depends on scripts 00-14. The purpose of this script is to compare genes
with differential gene expression, differential transcript usage, and
differential transcript expression. Run back in github docker. It is not
currently finished.

***Script 16: de_analysis/16_dataset_overview_figure.Rmd*** - The
purpose of this script is to provide a dataset overview to serve as
manuscript figure 1. It is dependent on scripts 00 and 01 to get sample
metadata. It takes less than 3 minutes to run and is not finished but is
mainly on the shelf for now.

#### Other Scripts

***dtu_analysis/functions.R*** - This is a function script for all Rmd
files in the dtu directory.

***de_analysis/de_functions.R*** - This is a function script for all Rmd
files in the de directory.

***Script 17: worm_analysis/17_dtu_whole_worm.Rmd*** - The purpose of
this script is to examine differential transcript usage by sex in
short-read c elegans data. I am especially interested in comparing this
DTU to my own mouse data to look for evolutionary conserved genes with
sex specific splicing. The entire script takes less than 10 minutes to
run. Must be run in at least docker 1.3 to include c elegans annotation.
This script does not depend on any other code, with the exception of
functions.R, which can be sourced.

***Script 18: worm_analysis/18_compare_worm_mouse.Rmd*** - The purpose
of this script is to compare worm and mouse genes with DTU across sex to
see if anything is conserved. Nothing was conserved so it was not super
useful.

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
