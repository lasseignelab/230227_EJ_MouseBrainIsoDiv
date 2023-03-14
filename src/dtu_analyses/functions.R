# the purpose of this script is to include all functions used in the dtu analyses

##################### dtu_region_region script #####################

# this function is for filtering genes that are tissue x in condition 1 or 2
filter_genes_or <- function(x) {
  # subset genes of interest
  sig_genes_subset <- filter(
    sig_isoform_genes,
    condition_1 == x | condition_2 == x
  )
  sig_genes_subset <- unique(sig_genes_subset$gene_id)
  # name object
  name <- substr(x, 1, 4)
  assign(paste0("sig_isoform_genes_", name),
         sig_genes_subset,
         envir = .GlobalEnv
  )
}

# filtering genes for each comparison, specifying tissue x for condition 1 and tissue y in condition 2
filter_genes_comp <- function(x, y) {
  # subset genes of interest
  sig_genes_subset <- filter(
    sig_isoform_genes,
    condition_1 == x & condition_2 == y
  )
  sig_genes_subset <- unique(sig_genes_subset$gene_id)
  # name object
  name_1 <- substr(x, 1, 4)
  name_2 <- substr(y, 1, 4)
  assign(paste0("sig_isoform_genes_", name_1, "_", name_2),
         sig_genes_subset,
         envir = .GlobalEnv
  )
}

# create volcano plot function, specifying tissue x for condition 1 and tissue y in condition 2
create_volcano_plot_comp <- function(x, y) {
  # create subset
  analyzed_subset <- dplyr::filter(
    cpm_switchlist_dexseq$isoformFeatures,
    condition_1 == x & condition_2 == y
  )
  # plot
  volcano <- ggplot(
    data = analyzed_subset,
    aes(x = dIF, y = -log10(isoform_switch_q_value))
  ) +
    geom_point(aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
               size = 2,
               alpha = 0.5,
               stroke = NA
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "magenta",
      linewidth = .7
    ) +
    geom_vline(
      xintercept = c(-0.1, 0.1),
      linetype = "dashed",
      color = "turquoise3",
      linewidth = .7
    ) +
    scale_color_manual("Signficant\nIsoform Switch",
                       labels = c("not significant", "significant"),
                       values = c("gray40", "limegreen")
    ) +
    labs(
      x = "dIF (Differential Isoform Fraction)",
      y = "-log10 (Isoform Switch q Value)"
    ) +
    theme_light() +
    theme(
      legend.text = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    ggtitle(paste0(
      "cerebellum", " vs. ", "cortex",
      " differentially used isoforms"
    ))
  # save
  ggsave(
    paste0(
      here("results", "plots", "DEXSeq_volcano"),
      "/", x, "_", y, "_volcano.png"
    ),
    plot = volcano, width = 6, height = 4
  )
}

# create function
run_plot_gprofiler <- function(x) {
  # get names
  name <- deparse(substitute(x))
  name <- substr(name, nchar(name) - 9 + 1, nchar(name))
  # remove numbers after decimal point
  genes_remove_decimal <- str_extract(
    x, "ENSMUSG..........."
  )
  # run gprofiler2
  gostres <- gost(
    query = genes_remove_decimal,
    organism = "mmusculus", ordered_query = FALSE,
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.05, correction_method = "g_SCS",
    domain_scope = "custom", custom_bg = bg_genes,
    numeric_ns = "", sources = NULL, as_short_link = FALSE
  )
  # plot
  gostplot(gostres, capped = TRUE, interactive = FALSE)
  # save
  ggsave(paste0(
    here("results", "plots", "gprofiler2"),
    "/", name, "_gostres.png"
  ), width = 6, height = 4)
}

##################### dtu_region_others script #####################

# create function for making switchlists and running dexseq on tissue x
make_switchlist_run_dexseq <- function(x) {
  # create design
  temp_design <- data.frame(
    sampleID = sample_collection_metadata$sample_id,
    condition = sample_collection_metadata$tissue == x
  )
  
  temp_design["condition"][temp_design["condition"] == TRUE] <- x
  temp_design["condition"][temp_design["condition"] == FALSE] <- "other"
  # create switchlist
  temp_switchlist <- importRdata(
    isoformCountMatrix = merged_counts_iso,
    isoformRepExpression = cpm_iso,
    designMatrix = temp_design,
    isoformExonAnnoation = here(
      "data", "nextflow", "results", "bambu", "extended_annotations.gtf"
    ),
    isoformNtFasta = here("data", "gffread", "isoform_sequences.fa"),
    showProgress = FALSE
  )
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run DEXSeq
  switchlist_analyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = FALSE
  )
  # rename object
  assign(paste0(x, "_switchlist_analyzed"),
         switchlist_analyzed,
         envir = .GlobalEnv
  )
  # save object
  saveRDS(switchlist_analyzed, here(
    "data", "switchlist_objects",
    paste0(x, "_switchlist_analyzed.Rds")
  ))
}
# create function to get gene symbols for all genes of a tissue x
get_gene_symbols <- function(x) {
  # assign name
  assign("temp_switchlist", get(paste0(x, "_switchlist_analyzed")))
  # pull shorter gene IDs
  temp_switchlist[["isoformFeatures"]]$shorter <-
    str_extract(temp_switchlist[["isoformFeatures"]]$gene_id,
                pattern = "ENSMUSG..........."
    )
  # get gene symbols from annotation dbi
  temp_gene_symbols <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(temp_switchlist[["isoformFeatures"]]$shorter),
    columns = c("SYMBOL", "ENSEMBL", "GENETYPE"), keytype = "ENSEMBL"
  )
  # add symbols and biotypes to object
  temp_switchlist[["isoformFeatures"]] <-
    left_join(temp_switchlist[["isoformFeatures"]],
              temp_gene_symbols,
              by = c("shorter" = "ENSEMBL")
    )
  # add to correct columns
  temp_switchlist[["isoformFeatures"]]$gene_name <-
    temp_switchlist[["isoformFeatures"]]$SYMBOL
  temp_switchlist[["isoformFeatures"]]$gene_biotype <-
    temp_switchlist[["isoformFeatures"]]$GENETYPE
  # remove extra columns
  temp_switchlist[["isoformFeatures"]]$shorter <- NULL
  temp_switchlist[["isoformFeatures"]]$SYMBOL <- NULL
  temp_switchlist[["isoformFeatures"]]$GENETYPE <- NULL
  # rename object
  assign(paste0(x, "_switchlist_analyzed"),
         temp_switchlist,
         envir = .GlobalEnv
  )
  saveRDS(temp_switchlist, here(
    "data", "switchlist_objects",
    paste0(x, "_switchlist_analyzed.Rds")
  ))
}
# make function for extracting significant genes
get_sig_genes <- function(x) {
  # subset features
  sig_features <- dplyr::filter(
    x$isoformFeatures,
    abs(dIF) > 0.1 & isoform_switch_q_value < 0.05
  )
  # get name
  name <- deparse(substitute(x))
  name <- substr(name, 1, 4)
  # assign object
  assign(paste0(name, "_sig_features"),
         sig_features,
         envir = .GlobalEnv
  )
  # pull genes
  sig_genes <- unique(sig_features$gene_id)
  # assign object
  assign(paste0(name, "_sig_genes"),
         sig_genes,
         envir = .GlobalEnv
  )
}

# create function
create_volcano_plot_region <- function(x) {
  # get name
  name <- deparse(substitute(x))
  name <- substr(name, 1, 4)
  # plot
  volcano_plot <- ggplot(
    data = x$isoformFeatures,
    aes(x = dIF, y = -log10(isoform_switch_q_value))
  ) +
    geom_point(aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
               size = 2,
               alpha = 0.5,
               stroke = NA
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "magenta",
      linewidth = .7
    ) +
    geom_vline(
      xintercept = c(-0.1, 0.1),
      linetype = "dashed",
      color = "turquoise3",
      linewidth = .7
    ) +
    scale_color_manual("Signficant\nIsoform Switch",
                       labels = c("not significant", "significant"),
                       values = c("gray40", "limegreen")
    ) +
    labs(
      x = "dIF (Differential Isoform Fraction)",
      y = "-log10 (Isoform Switch q Value)"
    ) +
    theme_light() +
    theme(
      legend.text = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    ggtitle(paste0(name, " differentially used isoforms"))
  # save
  ggsave(
    paste0(
      here("results", "plots", "DEXSeq_volcano"),
      "/", name, "_all_volcano.png"
    ),
    volcano_plot,
    width = 6, height = 4
  )
}

##################### dtu_region_sex script #####################

# make function to split out brain regions
split_region <- function(x) {
  # subset counts
  subset_counts <- subset(
    merged_counts_iso,
    select = c(
      "isoform_id",
      sample_collection_metadata$sample_id[
        sample_collection_metadata$tissue == x
      ]
    )
  )
  # subset cpm
  subset_cpm <- subset(
    cpm_iso,
    select = c(
      "isoform_id",
      sample_collection_metadata$sample_id[
        sample_collection_metadata$tissue == x
      ]
    )
  )
  # drop rows with 0 cpm
  subset_cpm <- subset_cpm[rowSums(subset_cpm[, -1]) != 0, ]
  # get name
  name <- substr(x, 1, 4)
  # assign objects
  assign(paste0(name, "_counts"), subset_counts, envir = .GlobalEnv)
  assign(paste0(name, "_cpm"), subset_cpm, envir = .GlobalEnv)
}

# make function for making a switchlist
make_switchlist <- function(x) {
  # get name
  name <- substr(x, 1, 4)
  # get objects for function
  assign("name_counts", get(paste0(name, "_counts")))
  assign("name_cpm", get(paste0(name, "_cpm")))
  # make design
  temp_design <- data.frame(
    sampleID = sample_collection_metadata$sample_id[
      sample_collection_metadata$tissue == x
    ],
    condition = sample_collection_metadata$sex[
      sample_collection_metadata$tissue == x
    ]
  )
  # make switchlist
  temp_switchlist <- importRdata(
    isoformCountMatrix = name_counts,
    isoformRepExpression = name_cpm,
    designMatrix = temp_design,
    isoformExonAnnoation = here(
      "data", "nextflow", "results", "bambu", "extended_annotations.gtf"
    ),
    isoformNtFasta = here("data", "gffread", "isoform_sequences.fa"),
    showProgress = FALSE
  )
  # assign variable
  assign(paste0(name, "_sex_switchlist"), temp_switchlist, envir = .GlobalEnv)
}

# make function to filter and run DEXSeq
filter_run_dexseq <- function(x) {
  # get name
  name <- substr(x, 1, 4)
  # assign name
  assign("temp_switchlist", get(paste0(name, "_sex_switchlist")))
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run DEXSeq
  switchlist_analyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = TRUE
  )
  # rename object
  assign(paste0(name, "_sex_switchlist_analyzed"),
         switchlist_analyzed,
         envir = .GlobalEnv
  )
  # save object
  saveRDS(switchlist_analyzed, here("data", "switchlist_objects",
                                    paste0(x, "_sex_switchlist_analyzed.Rds"))
  )  
}
# create function for sex split volcano plot for each brain region x
create_sex_volcano_plot <- function(x) {
  # get name
  name <- deparse(substitute(x))
  name <- substr(name, 1, 4)
  # plot
  volcano_plot <- ggplot(
    data = x$isoformFeatures,
    aes(x = dIF, y = -log10(isoform_switch_q_value))
  ) +
    geom_point(aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
               size = 2,
               alpha = 0.5,
               stroke = NA
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "magenta",
      linewidth = .7
    ) +
    geom_vline(
      xintercept = c(-0.1, 0.1),
      linetype = "dashed",
      color = "turquoise3",
      linewidth = .7
    ) +
    scale_color_manual("Signficant\nIsoform Switch",
                       labels = c("not significant", "significant"),
                       values = c("gray40", "limegreen")
    ) +
    labs(
      x = "dIF (Differential Isoform Fraction)",
      y = "-log10 (Isoform Switch q Value)"
    ) +
    theme_light() +
    theme(
      legend.text = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    ggtitle(paste0(name, " sex differentially used isoforms"))
  # save
  ggsave(
    paste0(
      here("results", "plots", "DEXSeq_volcano"),
      "/", name, "_sex_volcano.png"
    ),
    volcano_plot,
    width = 6, height = 4
  )
}
# make different function without reducing argument
filter_run_dexseq_noreduce <- function(x) {
  # get name
  name <- substr(x, 1, 4)
  # assign name
  assign("temp_switchlist", get(paste0(name, "_sex_switchlist")))
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run DEXSeq
  switchlist_analyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = FALSE
  )
  # rename object
  assign(paste0(name, "_sex_switchlist_analyzed"),
         switchlist_analyzed,
         envir = .GlobalEnv
  )
  # save object
  saveRDS(switchlist_analyzed, here("data", "switchlist_objects",
                                    paste0(x, "_sex_switchlist_analyzed.Rds"))
  )  
}

# make volcano plot function, but remove NAs
create_sex_volcano_plot_rm <- function(x) {
  # get name
  name <- deparse(substitute(x))
  name <- substr(name, 1, 4)
  # plot
  volcano_plot <- ggplot(
    data = x$isoformFeatures,
    aes(x = dIF, y = -log10(isoform_switch_q_value))
  ) +
    geom_point(aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
               size = 2,
               alpha = 0.5,
               stroke = NA
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "magenta",
      linewidth = .7
    ) +
    geom_vline(
      xintercept = c(-0.1, 0.1),
      linetype = "dashed",
      color = "turquoise3",
      linewidth = .7
    ) +
    scale_color_manual("Signficant\nIsoform Switch",
                       labels = c("not significant", "significant"),
                       values = c("gray40", "limegreen"),
                       na.translate = FALSE
    ) +
    labs(
      x = "dIF (Differential Isoform Fraction)",
      y = "-log10 (Isoform Switch q Value)"
    ) +
    theme_light() +
    theme(
      legend.text = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    ggtitle(paste0(name, " sex differentially used isoforms"))
  # save
  ggsave(
    paste0(
      here("results", "plots", "DEXSeq_volcano"),
      "/", name, "_sex_volcano.png"
    ),
    volcano_plot,
    width = 6, height = 4
  )
}

# write function to get significant isoforms
get_sig_isoforms <- function(x) {
  # get name of object
  name <- substr(x, 1, 4)
  # assign name
  assign(
    "temp_sex_switchlist_analyzed",
    get(paste0(name, "_sex_switchlist_analyzed"))
  )
  # get unique isoforms for that gene
  temp_sex_sig_isoforms <- unique(
    dplyr::filter(
      temp_sex_switchlist_analyzed$isoformFeatures,
      abs(dIF) > 0.1 & isoform_switch_q_value < 0.05
    )$isoform_id
  )
  # rename object
  assign(paste0(name, "_sex_sig_isoforms"),
         temp_sex_sig_isoforms,
         envir = .GlobalEnv
  )
}

# make simple function to get significant genes and cutting off decimals
get_cut_sig_genes <- function(x) {
  # get name
  name <- substr(x, 1, 4)
  # assign name
  assign(
    "temp_sex_switchlist_analyzed",
    get(paste0(name, "_sex_switchlist_analyzed"))
  )
  # subset features
  temp_sex_sig_features <- dplyr::filter(
    temp_sex_switchlist_analyzed$isoformFeatures,
    abs(dIF) > 0.1 & isoform_switch_q_value < 0.05
  )
  # assign object
  assign(paste0(name, "_sig_features"),
         temp_sex_sig_features,
         envir = .GlobalEnv
  )
  # pull genes
  sig_genes <- unique(temp_sex_sig_features$gene_id)
  sig_genes <- str_extract(sig_genes, "ENSMUSG...........")
  # assign object
  assign(paste0(name, "_sex_sig_genes"),
         sig_genes,
         envir = .GlobalEnv
  )
}
# make go analysis function that skips plotting if no result returned
run_gprofiler <- function(x) {
  # get name
  name <- substr(x, 1, 4)
  # assign name
  assign(
    "temp_sex_sig_genes",
    get(paste0(name, "_sex_sig_genes"))
  )
  # run gpfrofiler
  temp_sex_gostres <- gost(
    query = temp_sex_sig_genes, organism = "mmusculus", ordered_query = FALSE,
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.05, correction_method = "g_SCS", domain_scope = "custom",
    custom_bg = str_extract(
      unique(cere_sex_switchlist_analyzed$isoformFeatures$gene_id),
      pattern = "ENSMUSG..........."
    ),
    numeric_ns = "", sources = NULL, as_short_link = FALSE
  )
  # only save and plot if there are significant results
  if (!is.null(temp_sex_gostres)) {
    # save object
    assign(paste0(name, "_sex_gostres"),
           temp_sex_gostres,
           envir = .GlobalEnv
    )
    # plot object
    gostplot(temp_sex_gostres, capped = TRUE, interactive = FALSE)
  }
}
# write function to get symbols for each list
get_gene_symbols_sex <- function(x) {
  # get name of object
  name <- substr(x, 1, 4)
  # assign name
  assign("temp_sex_sig_genes", get(paste0(name, "_sex_sig_genes")))
  # get gene symbols from annotation dbi
  temp_sex_sig_gene_symbols <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = temp_sex_sig_genes,
    columns = c("SYMBOL", "GENENAME", "ENSEMBL"), keytype = "ENSEMBL"
  )
  # rename object
  assign(paste0(name, "_sex_sig_gene_symbols"),
         temp_sex_sig_gene_symbols,
         envir = .GlobalEnv
  )
}