# this script is for selecting out 5 or 3 samples of the same sex for a given
# and seeing if DTU genes detected changes much

library(tidyverse)
library(IsoformSwitchAnalyzeR)
library(ComplexHeatmap)
library(viridis)
library(org.Mm.eg.db)
library(gprofiler2)
library(styler)
library(lintr)
library(here)

source("src/dtu_analyses/calculate_cpm.R", local = knitr::knit_global())
source("src/dtu_analyses/functions.R", local = knitr::knit_global())

# pull out only female cortex and cerebellum
cort_cere_females <- sample_collection_metadata %>% filter(
  sex == "F" & (tissue == "cerebellum" | tissue == "cortex")
  )

# get only 3 mice for each
cort_cere_females_reduced <- cort_cere_females %>% filter(
  mouse_id == "mouse10" | mouse_id == "mouse11" | mouse_id == "mouse12")

# split counts - full
# subset counts
full_subset_counts <- subset(
  merged_counts_iso,
  select = c(
    "isoform_id", cort_cere_females$sample_id
  )
)
# subset cpm
full_subset_cpm <- subset(
  cpm_iso,
  select = c(
    "isoform_id", cort_cere_females$sample_id
  )
)
# drop rows with 0 cpm
full_subset_cpm <- full_subset_cpm[rowSums(full_subset_cpm[, -1]) != 0, ]

# make design - full
full_design <- data.frame(
    sampleID = cort_cere_females$sample_id,
    condition = cort_cere_females$tissue
  )
# make switchlist - full
full_switchlist <- importRdata(
    isoformCountMatrix = full_subset_counts,
    isoformRepExpression = full_subset_cpm,
    designMatrix = full_design,
    isoformExonAnnoation = here(
      "data", "nextflow", "results", "bambu", "extended_annotations.gtf"
    ),
    isoformNtFasta = here("data", "gffread", "isoform_sequences.fa"),
    showProgress = FALSE
  )



# repeat
# split counts - reduced
# subset counts
reduced_subset_counts <- subset(
  merged_counts_iso,
  select = c(
    "isoform_id", cort_cere_females_reduced$sample_id
  )
)
# subset cpm
reduced_subset_cpm <- subset(
  cpm_iso,
  select = c(
    "isoform_id", cort_cere_females_reduced$sample_id
  )
)
# drop rows with 0 cpm
reduced_subset_cpm <- reduced_subset_cpm[rowSums(reduced_subset_cpm[, -1]) != 0, ]

# make design - reduced
reduced_design <- data.frame(
  sampleID = cort_cere_females_reduced$sample_id,
  condition = cort_cere_females_reduced$tissue
)
# make switchlist - reduced
reduced_switchlist <- importRdata(
  isoformCountMatrix = reduced_subset_counts,
  isoformRepExpression = reduced_subset_cpm,
  designMatrix = reduced_design,
  isoformExonAnnoation = here(
    "data", "nextflow", "results", "bambu", "extended_annotations.gtf"
  ),
  isoformNtFasta = here("data", "gffread", "isoform_sequences.fa"),
  showProgress = FALSE
)

##### now run dexseq on both

# filter switchlist
full_switchlist <- preFilter(full_switchlist, geneExpressionCutoff = NULL)
# run DEXSeq
full_switchlist_analyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = full_switchlist,
  reduceToSwitchingGenes = TRUE)

# filter switchlist
reduced_switchlist <- preFilter(reduced_switchlist, geneExpressionCutoff = NULL)
# run DEXSeq
reduced_switchlist_analyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = reduced_switchlist,
  reduceToSwitchingGenes = TRUE)
  

# compare numbers

length(unique(full_switchlist_analyzed[["isoformFeatures"]]$gene_id))

length(unique(reduced_switchlist_analyzed[["isoformFeatures"]]$gene_id))

# having 3 samples does seem to show less features
