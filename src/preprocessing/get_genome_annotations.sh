#!/bin/bash

wget -P ../../data/gencode_annotations https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gtf.gz
wget -P ../../data/gencode_annotations https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz

gunzip ../../data/gencode_annotations/gencode.vM31.primary_assembly.annotation.gtf.gz
gunzip ../../data/gencode_annotations/GRCm39.primary_assembly.genome.fa.gz