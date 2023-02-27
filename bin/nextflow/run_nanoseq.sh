#!/bin/bash

module load Nextflow
module load Singularity

nextflow run nf-core/nanoseq -r 2.0.1 \
    --input samplesheet.csv \
    --protocol cDNA \
    --flowcell FLO-MIN106 \
    --kit SQK-PCB109 \
    --skip_basecalling \
    --skip_demultiplexing \
    --skip_differential_analysis \
    -profile cheaha \
    -c custom.config