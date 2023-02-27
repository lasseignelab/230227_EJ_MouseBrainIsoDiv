#code for reading and preparing data for environment

#load in packages
library(tidyverse)

#read in metadata
sample_collection_metadata <- read.csv("../../doc/sample_collection_metadata.csv")

#read in counts data
merged_counts <- read.table("../../data/nextflow/results/bambu/counts_transcript.txt", header = TRUE)
merged_counts_iso <- merged_counts[,-2]
merged_counts_noid <- merged_counts[,-c(1,2)]

colnames(merged_counts_iso) <- c("isoform_id", sample_collection_metadata$sample_id)
colnames(merged_counts_noid) <- c(sample_collection_metadata$sample_id)

#get cpm
cpm <- do.call(cbind, lapply(1:ncol(merged_counts_noid), function(i) {
  merged_counts_noid[i]*1e6 / sum(merged_counts_noid[i])
}))
cpm_iso <- cbind(merged_counts[,1], cpm)
colnames(cpm_iso) <- colnames(merged_counts_iso)
