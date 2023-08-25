library(tidyverse)
library(biomaRt)

# read in counts data
gene_exp_counts <- read.table("data/nextflow/bambu/counts_gene.txt",
                              header = TRUE)

# calculate cpm
gene_exp_cpm <-
  do.call(cbind, lapply(seq_len(ncol(gene_exp_counts)), function(i) {
    gene_exp_counts[i] * 1e6 / sum(gene_exp_counts[i])
  }))

# remove novel genes
gene_exp_cpm <- gene_exp_cpm[!grepl("gene", rownames(gene_exp_cpm)),]

# pull IDs
ensembl_ids <- str_extract(rownames(gene_exp_cpm), "ENSMUSG...........")

# rename rownames
rownames(gene_exp_cpm) <- ensembl_ids

# pull rows with symbols
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

matched_symbols <- getBM(filters = "ensembl_gene_id",
              attributes = c("ensembl_gene_id", "mgi_symbol"), values = ensembl_ids, 
              mart = mart)

duplicated_symbols <- matched_symbols$mgi_symbol[duplicated(matched_symbols$mgi_symbol)]

duplicated_symbols <- unique(duplicated_symbols)[-1] # remove all the blanks

matched_symbols <- matched_symbols[!duplicated(matched_symbols$mgi_symbol),] # remove duplicates, need to use ensembl ID for those

ensembl_order <- gene_exp_cpm[matched_symbols$ensembl_gene_id,]

symbols_cpm <- ensembl_order

rownames(symbols_cpm) <- matched_symbols$mgi_symbol

combined_cpm <- rbind(gene_exp_cpm, symbols_cpm)

colnames(combined_cpm) <- substr(colnames(combined_cpm), 1, 8)

# save this object to shiny app directory
saveRDS(combined_cpm, "src/shiny_app/data/combined_cpm.Rds")

# read in sample collection metadata
sample_collection_metadata <- read_csv("doc/sample_collection_metadata.csv")

# reorder metadata
order <- colnames(combined_cpm)
order <- substr(order, 1, 8)

sample_collection_metadata_reorder <- 
  match(order,sample_collection_metadata$sample_id)
sample_collection_metadata <- 
  sample_collection_metadata[sample_collection_metadata_reorder,]

# save this object to shiny app directory
saveRDS(sample_collection_metadata, "src/shiny_app/data/sample_collection_metadata.Rds")
