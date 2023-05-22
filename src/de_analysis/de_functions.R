# function for getting results from a produced DESeq2 object based on tissues
get_de_results <- function(tissue1, tissue2) {
  temp_results <- results(dds,
                          contrast = c("tissue", tissue1, tissue2),
                          independentFiltering = TRUE, alpha = 0.05,
                          pAdjustMethod = "BH", parallel = TRUE
  )
  saveRDS(temp_results, here("data", "deseq2_data", paste0(
    tissue1, "_",
    tissue2,
    "_results.Rds"
  )))
  assign(paste0(tissue1, "_", tissue2, "_res"), temp_results, env = .GlobalEnv)
  
  sig_temp_results <- subset(temp_results, padj < 0.05)
  
  sig_genes <- rownames(sig_temp_results[
    order(sig_temp_results$log2FoldChange),
  ])
  
  print(paste0(
    tissue1, " versus ", tissue2, " had ",
    length(sig_genes), " significant DE genes"
  ))
  
  
  write.table(sig_genes,
              file = here(
                "results", "de_genes",
                paste0(tissue1, "_", tissue2, ".csv")
              ),
              row.names = FALSE, quote = FALSE, col.names = FALSE
  )
}