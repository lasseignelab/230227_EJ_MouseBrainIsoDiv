Nanopore long-read RNA sequencing shows region-specific sex differences
in wild-type mouse brain mRNA isoform expression and usage
================
2023-12-01

## Authors

**Emma F. Jones, Timothy C. Howton, Victoria L. Flanary, Amanda D.
Clark, and Brittany N. Lasseigne.**

The University of Alabama at Birmingham, Heersink School of Medicine,
Department of Cell, Developmental and Integrative Biology

## Project Overview

The purpose of this project is to analyze Oxford Nanopore RNA sequencing
from four wild-type mouse brain regions balanced for sex, to assay
isoform usage and expression differences across brain region and sex.

![We extracted RNA from 40 samples, four brain regions each from 10 mice
balanced for sex with and sequenced on an ONT GridION. We processed this
data with the nf-core nanoseq pipeline and used DESeq2,
IsoformSwitchAnalyzeR, and custom scripts for downstream
analyses](https://github.com/lasseignelab/230227_EJ_MouseBrainIsoDiv/blob/main/src/shiny_app/www/graphical_abstract.png)
Preprint DOI  
Docker DOI  
Zenodo DOI  
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE246705>

You can also visit our Shiny Application that accompanies this analysis:

<https://lasseignelab.shinyapps.io/mouse_brain_iso_div/>

## Scripts

#### Preprocessing

    ## src/preprocessing
    ## ├── 00_run_nanoseq.sh
    ## ├── 02_get_genome_annotations.sh
    ## ├── 03_create_isoform_fa.sh
    ## ├── 16_dataset_overview_figure1.Rmd
    ## ├── 20010631.c0158.err.txt
    ## ├── 20010631.c0158.out.txt
    ## └── README

***Script 00: preprocessing/00_run_nanoseq.sh*** - This script can be
skipped if you are working from counts files. Does not require a docker.
You need raw input fasta files, a genome file, and genome annotation
file. I used GENCODE mouse release M31.

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

***Script 16: preprocessing/16_dataset_overview_figure.Rmd*** - The
purpose of this script is to provide a dataset overview to serve as
manuscript figure 1. It is dependent on scripts 00 and 01 to get sample
metadata. Must run in github docker 1.7 or above.

#### DTU Analyses

    ## src/dtu_analyses
    ## ├── 01_calculate_cpm.R
    ## ├── 04_pca_eda.Rmd
    ## ├── 05_dtu_region_region.Rmd
    ## ├── 06_dtu_region_others.Rmd
    ## ├── 07_dtu_region_sex.Rmd
    ## ├── 08_dtu_neuro_diseases.Rmd
    ## ├── 09_dtu_isoform_switching.Rmd
    ## ├── 14_protein_domain_info.Rmd
    ## ├── README
    ## ├── functions.R
    ## └── size_power.R

***dtu_analysis/functions.R*** - This is a function script for all Rmd
files in the dtu directory.

***Script 01: dtu_analyses/01_calculate_cpm.R*** - This script is very
short and depends on having transcript counts available either by
downloading or running the nanoseq pipeline. Takes less than one minute.

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

***Script 14: dtu_analyses/14_protein_domain_info.Rmd*** - This script
depends on scripts 00-13, and also includes a shell script for running
perl code. The purpose of this script is to extract nucleotide and amino
acid sequences and run pfam to annotate the protein domains. NEEDS TO
RUN IN PFAM DOCKER!

#### DGE and DTE Analyses

    ## src/de_analyses
    ## ├── 10_DESeq2_region_region.Rmd
    ## ├── 11_DESeq2_region_others.Rmd
    ## ├── 12_DESeq2_region_sex.Rmd
    ## ├── 13_incorporate_de_results.Rmd
    ## ├── 15_compare_DTU_DGE.Rmd
    ## ├── 17_compare_region_pairwise_results.Rmd
    ## ├── 18_enrichment_analysis.Rmd
    ## ├── 19_create_fig_2.Rmd
    ## ├── README
    ## ├── de_functions.R
    ## └── stacked_barplot_functions.R

***de_analyses/de_functions.R*** - This is a function script for all Rmd
files in the de directory.

***de_analyses/stacked_barplot_functions.R*** - This is a function
script for specifically building stacked barplots for figure 2.

***Script 10: de_analyses/10_DESeq2_region_region.Rmd*** - This script
is dependent on script 00/01 or having gene and transcript level count
data with metadata in your /data/ directory and the de functions script.
The purpose of this script is to run DESeq2 across brain regions at the
gene and transcript level. It takes less than 5 minutes to run.

***Script 11: de_analyses/11_DESeq2_region_others.Rmd*** - This script
is dependent on script 00/01 or having gene and transcript level count
data with metadata in your /data/ directory and the de functions script.
The purpose of this script is to run DESeq2 for each brain region
compared to an aggregate of other brain regions at the gene and
transcript level. It takes about 10 minutes to run.

***Script 12: de_analyses/12_DESeq2_region_sex.Rmd*** - This script is
dependent on script 00/01 or having gene and transcript level count data
with metadata in your /data/ directory and the de functions script. The
purpose of this script is to compare expression at the gene and
transcript level across sexes within brain regions. It takes less than
10 minutes to run.

***Script 13: de_analyses/13_incorporate_de_results.Rmd*** - This script
is fully dependent on scripts 01-09 and the de_analysis DESeq2 scripts
10-12 and the de functions script. The purpose of this script is to wrap
in all DESeq2 significance values into my isoformSwitchAnalyzeR objects
and subsequently, plots. This enables us to plot them all with
significance values for DGE, DTE, and DTU. This script takes 12 minutes
to run.

***Script 15: de_analyses/15_compare_DTU_DGE.Rmd*** - This script
depends on scripts 00-14. The purpose of this script is to compare genes
with differential gene expression, differential transcript usage, and
differential transcript expression. Must run in docker 1.7 to have all
packages.

***Script 17: de_analyses/17_compare_region_pairwise_results.Rmd*** -
The purpose of this script is to determine how many DTU genes or DEGs
are specific to a study design. I used Jaccard Similarity to see how
similar the results are to each other and create supp figure 3. Must run
in github docker 1.7 or above.

***Script 18: de_analyses/18_enrichment_analysis.Rmd*** - The purpose of
this script is to perform functional enrichment analysis with the
package gprofiler2 on the DGE and DTE genes. This script is fully
dependent on scripts 00-12. Must run in github docker 1.7 or above.

***Script 19: de_analyses/19_create_fig_2.Rmd*** - The purpose of this
script is to create supplementary figure 2. This script is fully
dependent on scripts 00-12. Must run in github docker 1.7 or above.

#### Shiny Application scripts

    ## src/shiny_support
    ## ├── README
    ## └── calculate_cpm_convert_ids.R

***shiny_support/calculate_cpm_convert_ids.R*** - This script enables
the gene IDs for the Shiny App to be either ENSEMBL IDs or gene symbols
and does some other needed reformatting of the data for the application
to work.

    ## src/shiny_app
    ## ├── README
    ## ├── app.R
    ## ├── complete_switchlists
    ## │   ├── region_all_list_orf_de_pfam.Rds
    ## │   ├── region_region_orf_de_pfam.Rds
    ## │   └── region_sex_list_orf_de_pfam.Rds
    ## ├── server.R
    ## ├── ui.R
    ## └── www
    ##     ├── DTU_example.png
    ##     ├── code.js
    ##     ├── favicon.ico
    ##     ├── graphical_abstract.png
    ##     ├── logo_only.png
    ##     └── style.css

***shiny_app/app.R*** - This script runs the entire Shiny application.

***shiny_app/server.R*** - The purpose of this script is to run the
back-end code needed for generating plots and data for the Shiny
application.

***shiny_app/ui.R*** - The purpose of this script is to run the user
interface for shiny web application.

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

## Funding

NHGRI R00HG009678 (PI: Lasseigne)  
Pittman Scholar (PI: Lasseigne)  
UAB Lab Startup funds (PI: Lasseigne)

## Acknowledgements

We acknowledge all current and past members of the Lasseigne Lab for
their thoughtful feedback, especially Tabea M. Soelter, Jordan H.
Whitlock, Vishal H. Oza, and Elizabeth J. Wilk. We would like to thank
the UAB Biological Data Sciences (UAB-BDS) core for its expertise,
institutional support, and maintenance of the nf-core nanoseq pipeline
and docker/singularity container documentation.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.

[![MIT
License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)
