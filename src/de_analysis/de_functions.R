# function for running deseq2 across brain regions looking at sex
deseq2_all_regions_sex <-
  function(counts_table,
           metadata,
           object_save_path,
           results_save_path,
           level) {
    # create deseq data set
    dds_region <- DESeqDataSetFromMatrix(
      countData = round(counts_table),
      colData = metadata,
      design = ~sex
    )
    # filter genes
    keep <- rowSums(counts(dds_region)) > 1
    dds_region <- dds_region[keep, ]
    nrow(dds_region)
    # run deseq2
    dds_region <- estimateSizeFactors(dds_region)
    dds_region <- estimateDispersions(dds_region)
    dds_region <- nbinomWaldTest(dds_region, maxit = 5000)
    # save object to global environment
    assign("sex_dds", dds_region, envir = .GlobalEnv)
    # pull results
    temp_results <- results(
      dds_region,
      independentFiltering = TRUE,
      alpha = 0.05,
      pAdjustMethod = "BH",
      parallel = TRUE
    )
    # save results output as well
    assign("sex_res", temp_results, env = .GlobalEnv)
    # save results object
    saveRDS(
      temp_results,
      paste0(
        object_save_path,
        "/all_regions_sex_",
        level,
        "_results.Rds"
      )
    )
    # filter for genes with an adj pval of less than 0.05 and log2FC over 1.5
    sig_temp_results <- subset(temp_results, padj < 0.05)
    sig_temp_results <-
      subset(sig_temp_results, abs(log2FoldChange) > 1.5)

    sig_genes <- rownames(
      sig_temp_results[order(sig_temp_results$log2FoldChange), ]
    )
    # print number of significant genes
    print(paste0(
      "M vs F had ",
      length(sig_genes),
      " significant DE ",
      level,
      "s by sex"
    ))
    # save csv of significant genes
    write.table(
      sig_genes,
      file = paste0(results_save_path, "/all_regions_sex", ".csv"),
      row.names = FALSE,
      quote = FALSE,
      col.names = FALSE
    )
  }

# function for running deseq2 on a single brain region
deseq2_single_region_sex <-
  function(counts_table,
           metadata,
           region,
           object_save_path,
           results_save_path,
           level) {
    # split out brain region from counts
    save_ids <- metadata$sample_id[metadata$tissue == region]
    # create deseq data set
    dds_region <- DESeqDataSetFromMatrix(
      countData = round(counts_table[, save_ids]),
      colData = metadata[metadata$tissue == region, ],
      design = ~sex
    )
    # filter genes
    keep <- rowSums(counts(dds_region)) > 1
    dds_region <- dds_region[keep, ]
    nrow(dds_region)
    # run deseq2
    dds_region <- estimateSizeFactors(dds_region)
    dds_region <- estimateDispersions(dds_region)
    dds_region <- nbinomWaldTest(dds_region, maxit = 5000)
    # save object to global environment
    assign(paste0("dds_", region), dds_region, envir = .GlobalEnv)
    # pull results
    temp_results <- results(
      dds_region,
      independentFiltering = TRUE,
      alpha = 0.05,
      pAdjustMethod = "BH",
      parallel = TRUE
    )
    # save results output as well
    assign(paste0(region, "_sex_res"), temp_results, env = .GlobalEnv)
    # save results object
    saveRDS(
      temp_results,
      paste0(
        object_save_path, "/",
        region, "_sex_", level,
        "_results.Rds"
      )
    )
    # filter for genes with an adj pval of less than 0.05 and log2FC over 1.5
    sig_temp_results <- subset(temp_results, padj < 0.05)
    sig_temp_results <-
      subset(sig_temp_results, abs(log2FoldChange) > 1.5)

    sig_genes <- rownames(
      sig_temp_results[order(sig_temp_results$log2FoldChange), ]
    )
    # print number of significant genes
    print(paste0(
      region,
      " had ",
      length(sig_genes),
      " significant DE ",
      level,
      "s by sex"
    ))
    # save csv of significant genes
    write.table(
      sig_genes,
      file = paste0(results_save_path, "/", region, "_sex", ".csv"),
      row.names = FALSE,
      quote = FALSE,
      col.names = FALSE
    )
  }

# function for getting results from a produced DESeq2 object based on tissues
get_de_results <-
  function(tissue1,
           tissue2,
           data_save_path,
           results_save_path) {
    temp_results <- results(
      dds,
      contrast = c("tissue", tissue1, tissue2),
      independentFiltering = TRUE,
      alpha = 0.05,
      pAdjustMethod = "BH",
      parallel = TRUE
    )
    saveRDS(
      temp_results,
      paste0(
        data_save_path, "/",
        tissue1, "_",
        tissue2,
        "_results.Rds"
      )
    )
    assign(paste0(tissue1, "_", tissue2, "_res"),
      temp_results,
      env = .GlobalEnv
    )

    # filter for signficant genes with a log2FC cutoff of 1.5
    sig_temp_results <- subset(temp_results, padj < 0.05)
    sig_temp_results <-
      subset(sig_temp_results, abs(log2FoldChange) > 1.5)

    sig_genes <- rownames(
      sig_temp_results[order(sig_temp_results$log2FoldChange), ]
    )

    print(paste0(
      tissue1,
      " versus ",
      tissue2,
      " had ",
      length(sig_genes),
      " significant DE genes"
    ))

    write.table(
      sig_genes,
      file = paste0(results_save_path, "/", tissue1, "_", tissue2, ".csv"),
      row.names = FALSE,
      quote = FALSE,
      col.names = FALSE
    )
  }

# function for getting same deseq results but from DTE and not DGE
get_de_results_transcripts <-
  function(tissue1,
           tissue2,
           data_save_path,
           results_save_path) {
    temp_results <- results(
      dds_transcripts,
      contrast = c("tissue", tissue1, tissue2),
      independentFiltering = TRUE,
      alpha = 0.05,
      pAdjustMethod = "BH",
      parallel = TRUE
    )
    saveRDS(
      temp_results,
      paste0(
        data_save_path,
        "/",
        tissue1,
        "_",
        tissue2,
        "_transcripts_results.Rds"
      )
    )
    assign(paste0(tissue1, "_", tissue2, "_transcript_res"),
      temp_results,
      env = .GlobalEnv
    )

    sig_temp_results <- subset(temp_results, padj < 0.05)
    sig_temp_results <-
      subset(sig_temp_results, abs(log2FoldChange) > 1.5)

    sig_transcripts <- rownames(
      sig_temp_results[order(sig_temp_results$log2FoldChange), ]
    )

    print(paste0(
      tissue1,
      " versus ",
      tissue2,
      " had ",
      length(sig_transcripts),
      " significant DE transcripts"
    ))

    write.table(
      sig_transcripts,
      file = paste0(results_save_path, "/", tissue1, "_", tissue2, ".csv"),
      row.names = FALSE,
      quote = FALSE,
      col.names = FALSE
    )
  }

# function for reformatting results so they match the switchlist data frame
format_deseq_results <- function(condition1, condition2) {
  assign(
    "condition1_condition2_res",
    get(paste0(condition1, "_", condition2, "_res"))
  )

  condition1_condition2_padj <- condition1_condition2_res[, 6]

  condition1_condition2_padj <- data.frame(
    "gene_id" = rownames(condition1_condition2_res),
    "padj" = condition1_condition2_padj
  )

  condition1_condition2_order <- order_needed %>%
    filter(condition_1 == condition1, condition_2 == condition2)

  condition1_condition2_both <-
    left_join(
      condition1_condition2_order,
      condition1_condition2_padj
    )

  final_gene_column <-
    c(final_gene_column, condition1_condition2_both$padj)

  assign("final_gene_column", final_gene_column, env = .GlobalEnv)
}

# function for reformatting results so they match the switchlist data frame, but
# at the transcript level instead of the gene level
format_deseq_results_dte <- function(condition1, condition2) {
  assign(
    "condition1_condition2_res",
    get(paste0(condition1, "_", condition2, "_dte_res"))
  )

  condition1_condition2_padj <- condition1_condition2_res[, 6]

  condition1_condition2_padj <- data.frame(
    "isoform_id" = rownames(condition1_condition2_res),
    "padj" = condition1_condition2_padj
  )

  condition1_condition2_order <- order_needed_transcripts %>%
    filter(condition_1 == condition1, condition_2 == condition2)

  condition1_condition2_both <-
    left_join(
      condition1_condition2_order,
      condition1_condition2_padj
    )

  final_transcript_column <-
    c(final_transcript_column, condition1_condition2_both$padj)

  assign("final_transcript_column", final_transcript_column, env = .GlobalEnv)
}

#######
# function for doing deseq2 on a single brain region
deseq2_single_region <- function(counts_table,
                                 metadata,
                                 region,
                                 object_save_path,
                                 results_save_path,
                                 level) {
  # add extra column
  metadata <- mutate(metadata, region = tissue == region)
  # create deseq data set
  dds_region <- DESeqDataSetFromMatrix(
    countData = round(counts_table),
    colData = metadata,
    design = ~region
  )
  # filter genes
  keep <- rowSums(counts(dds_region)) > 1
  dds_region <- dds_region[keep, ]
  nrow(dds_region)
  # run deseq2
  dds_region <- estimateSizeFactors(dds_region)
  dds_region <- estimateDispersions(dds_region)
  dds_region <- nbinomWaldTest(dds_region, maxit = 5000)
  # save object to global environment
  assign(paste0("dds_", region), dds_region, envir = .GlobalEnv)
  # pull results
  temp_results <- results(
    dds_region,
    independentFiltering = TRUE,
    alpha = 0.05,
    pAdjustMethod = "BH",
    parallel = TRUE
  )
  # save results output as well
  assign(paste0(region, "_res"), temp_results, env = .GlobalEnv)
  # save results object
  saveRDS(
    temp_results,
    paste0(
      object_save_path, "/",
      region, "_", level,
      "_results.Rds"
    )
  )
  # filter for genes with an adjusted pval of less than 0.05 and log2FC over 1.5
  sig_temp_results <- subset(temp_results, padj < 0.05)
  sig_temp_results <-
    subset(sig_temp_results, abs(log2FoldChange) > 1.5)

  sig_genes <- rownames(
    sig_temp_results[order(sig_temp_results$log2FoldChange), ]
  )
  # print number of significant genes
  print(paste0(
    region, " had ",
    length(sig_genes), " significant DE ", level, "s"
  ))
  # save csv of significant genes
  write.table(
    sig_genes,
    file = paste0(results_save_path, "/", region, ".csv"),
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE
  )
}

# function for incorporating both gene and transcript expression deseq results
# with a switchanalayzerlist
incorporate_deseq_results <-
  function(deseq2_gene_results,
           deseq2_tx_results,
           switchlist_obj) {
    # get gene results table
    gene_padj <- deseq2_gene_results[, 6]
    # pull just gene id and adjusted p values
    gene_padj <- data.frame(
      "gene_id" =
        rownames(deseq2_gene_results),
      "padj" = gene_padj
    )
    # reorder here and make longer
    gene_joined <-
      left_join(switchlist_obj[["isoformFeatures"]], gene_padj)
    # add the matched data
    switchlist_obj$isoformFeatures$gene_q_value <-
      gene_joined$padj
    # get tx results table
    transcript_padj <- deseq2_tx_results[, 6]
    # pull just gene id and adjusted p values
    transcript_padj <- data.frame(
      "isoform_id" =
        rownames(deseq2_tx_results),
      "padj" = transcript_padj
    )
    # reorder here and make longer
    tx_joined <-
      left_join(switchlist_obj[["isoformFeatures"]], transcript_padj)
    # add the matched data
    switchlist_obj$isoformFeatures$iso_q_value <- tx_joined$padj
    # return object
    return(switchlist_obj)
  }
