#!/bin/bash

mkdir -p ../../data/gffread
gffread -w ../../data/gffread/isoform_sequences.fa -g ../../data/gencode_annotations/GRCm39.primary_assembly.genome.fa ../../data/nextflow/bambu/extended_annotations.gtf 