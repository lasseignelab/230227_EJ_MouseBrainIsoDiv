#!/bin/bash

gffread -w ../../data/merged_fastq_annotations/isoform_sequences.fa -g ../../data/gencode_annotations/GRCm39.primary_assembly.genome.fa ../../data/merged_fastq_annotation/extended_annotations.gtf 