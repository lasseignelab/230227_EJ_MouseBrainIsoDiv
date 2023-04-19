# code for reading and preparing data for environment

# load in packages
library(tidyverse)
library(here)

# proc time
ptm <- proc.time()

# read in metadata
sample_collection_metadata <- read.csv(
  here("doc", "sample_collection_metadata.csv")
)

# read in counts data
merged_counts <- read.table(
  here("data", "nextflow", "bambu", "counts_transcript.txt"),
  header = TRUE
)
merged_counts_iso <- merged_counts[, -2]
merged_counts_noid <- merged_counts[, -c(1, 2)]

colnames(merged_counts_iso) <- c(
  "isoform_id",
  sample_collection_metadata$sample_id
)
colnames(merged_counts_noid) <- c(sample_collection_metadata$sample_id)

# get cpm
cpm <- do.call(cbind, lapply(seq_len(ncol(merged_counts_noid)), function(i) {
  merged_counts_noid[i] * 1e6 / sum(merged_counts_noid[i])
}))
cpm_iso <- cbind(merged_counts[, 1], cpm)
colnames(cpm_iso) <- colnames(merged_counts_iso)

# save RDS
save.image(file = here("data", "cpm_out", "cpm_counts_metadata.RData"))

# proc time
fptm <- proc.time() - ptm
fptm[3]/60

#session info
sessionInfo()

#R version 4.2.3 (2023-03-15)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.2 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] here_1.0.1      lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.1     purrr_1.0.1    
# [7] readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.2   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] rstudioapi_0.14  magrittr_2.0.3   hms_1.1.3        tidyselect_1.2.0 munsell_0.5.0   
# [6] timechange_0.2.0 colorspace_2.1-0 R6_2.5.1         rlang_1.1.0      fansi_1.0.4     
# [11] tools_4.2.3      grid_4.2.3       gtable_0.3.3     utf8_1.2.3       cli_3.6.1       
# [16] withr_2.5.0      rprojroot_2.0.3  lifecycle_1.0.3  tzdb_0.3.0       vctrs_0.6.1     
# [21] glue_1.6.2       stringi_1.7.12   compiler_4.2.3   pillar_1.9.0     generics_0.1.3  
# [26] scales_1.2.1     pkgconfig_2.0.3 