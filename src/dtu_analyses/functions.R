# Function script for use with the EJMouseBrainIsoDiv project and was written by
# Emma Jones. The purpose of this script is to include all functions used in DTU
# analyses. The input "tissue" or "region" is a character vector denoting brain
# region. This file is broken up into sections with functions for each
#individual numbered script.

######################### 04_pca_eda script ####################################

# this function is for PCA plotting
plot_pca <- function(metadata, firstPC, secondPC, color, shape) {
  color_len <- deparse(substitute(color))
  shape_len <- deparse(substitute(shape))
  color_class <- class(select(metadata, color_len)[[1]])
  p <- ggplot(metadata, aes(x = {{ firstPC }}, y = {{ secondPC }})) +
    geom_point(aes(color = {{ color }}, shape = {{ shape }}), size = 3) +
    theme_bw(base_size = 12) +
    labs(
      x = paste0(
        deparse(substitute(firstPC)),
        ": ", 
        round(var_explained[
          as.numeric(substr(deparse(substitute(firstPC)), 3, 4))
        ] * 100, 1), 
        "%"),
      y = paste0(
        deparse(substitute(secondPC)),
        ": ", 
        round(var_explained[
          as.numeric(substr(deparse(substitute(secondPC)), 3, 4))
        ] * 100, 1), 
        "%")
    ) 
    
  if (length(unique(metadata[[color_len]])) > 10
      && color_class == "character"
      | length(unique(metadata[[shape_len]])) > 10) {
    p +
      geom_text(
        data = metadata,
        label = metadata[[color_len]],
        nudge_x = 0.25, nudge_y = 0.25,
        check_overlap = TRUE
      ) +
      theme(legend.position = "none") 
  } else {  
    p +
      theme(legend.position = "top")
  }
}

######################### 05-07 scripts ##########################
# This function is for creating the swithlist.
make_switchlist_saturn <- function(isoformCountMatrix,
                                   isoformRepExpression,
                                   designMatrix,
                                   isoformExonAnnoation,
                                   isoformNtFasta,
                                   reduceToSwitchingGenes) {
  # create switchlist
  temp_switchlist <- importRdata(
    isoformCountMatrix = isoformCountMatrix,
    isoformRepExpression = isoformRepExpression,
    designMatrix = designMatrix,
    isoformExonAnnoation = isoformExonAnnoation,
    isoformNtFasta = isoformNtFasta,
    showProgress = FALSE
  )
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run satuRn
  switchlist_analyzed <- isoformSwitchTestSatuRn(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = reduceToSwitchingGenes
  )
  return(switchlist_analyzed)
}

# This function is for creating the swithlist with DEXSeq.
make_switchlist_dexseq <- function(isoformCountMatrix,
                                   isoformRepExpression,
                                   designMatrix,
                                   isoformExonAnnoation,
                                   isoformNtFasta,
                                   reduceToSwitchingGenes) {
  # create switchlist
  temp_switchlist <- importRdata(
    isoformCountMatrix = isoformCountMatrix,
    isoformRepExpression = isoformRepExpression,
    designMatrix = designMatrix,
    isoformExonAnnoation = isoformExonAnnoation,
    isoformNtFasta = isoformNtFasta,
    showProgress = FALSE
  )
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run satuRn
  switchlist_analyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = reduceToSwitchingGenes
  )
  return(switchlist_analyzed)
}

# this function is for getting gene symbols for all genes of a tissue
# save_path for this should be here("data", "switchlist_objects")



# Rewriting get_gene_symbols
get_gene_symbols <- function(switchlist_obj) {
  # pull shorter gene IDs
  switchlist_obj[["isoformFeatures"]]$shorter <-
    str_extract(switchlist_obj[["isoformFeatures"]]$gene_id,
                pattern = "ENSMUSG..........."
    )
  # get gene symbols from annotation dbi
  temp_gene_symbols <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(switchlist_obj[["isoformFeatures"]]$shorter),
    columns = c("SYMBOL", "ENSEMBL", "GENETYPE"), keytype = "ENSEMBL"
  )
  # add symbols and biotypes to object
  switchlist_obj[["isoformFeatures"]] <-
    left_join(switchlist_obj[["isoformFeatures"]],
              temp_gene_symbols,
              by = c("shorter" = "ENSEMBL")
    )
  # add to correct columns
  switchlist_obj[["isoformFeatures"]]$gene_name <-
    switchlist_obj[["isoformFeatures"]]$SYMBOL
  switchlist_obj[["isoformFeatures"]]$gene_biotype <-
    switchlist_obj[["isoformFeatures"]]$GENETYPE
  # remove extra columns
  switchlist_obj[["isoformFeatures"]]$shorter <- NULL
  switchlist_obj[["isoformFeatures"]]$SYMBOL <- NULL
  switchlist_obj[["isoformFeatures"]]$GENETYPE <- NULL
  return(switchlist_obj)
}

# this function is for filtering genes for each comparison
# specifying tissue1 for condition 1 and tissue2 in condition 2
filter_genes <- function(comparisons, switchlist_obj, sig_isoform_features) {
  sig_overlap_list <- list()
  for (comparison in comparisons) {
    # set up conditions for filtering
    con_1 <- str_split_1(comparison, "_")[1]
    con_2 <- str_split_1(comparison, "_")[2]
    # pull just genes
    sig_isoform_genes <- sig_isoform_features %>%
      distinct(gene_id, condition_1, condition_2, .keep_all = TRUE)
    # subset genes of interest
    sig_genes_subset <- filter(
      sig_isoform_genes,
      condition_1 == con_1 & condition_2 == con_2
    )
    sig_genes_subset <- unique(sig_genes_subset$gene_name)
    sig_genes_subset <- na.omit(sig_genes_subset)
    sig_overlap_list[[comparison]] <- sig_genes_subset
  }
  return(sig_overlap_list)
}
  
create_volcano_plot <- function(switchlist, condition1=NULL, condition2=NULL) {
  if(!is.null(condition1) & !is.null(condition2)) {
    # create subset
    switchlist_data <- dplyr::filter(
      switchlist$isoformFeatures,
      condition_1 == condition1 & condition_2 == condition2
    )
  } else {
    switchlist_data <- switchlist$isoformFeatures
  }
  # plot
  volcano <- ggplot(
    data = switchlist_data,
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
    ggtitle(paste0(
      condition1, " vs. ", condition2,
      " differentially used isoforms"
    ))
  volcano
}

# this function is for runnning gprofiler
run_plot_gprofiler <- function(gene_list,
                               name,
                               save_path) {
  # run gprofiler2
  gostres <- gost(
    query = gene_list,
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
    save_path,
    "/", name, "_gostres.png"
  ), width = 6, height = 4)
}

######################### 06_dtu_region_others script ##########################

# this function is for making switchlists and running saturn on a tissue
# must specify a save_path as well, use here() for easier path
make_switchlist_run_saturn <- function(tissue, save_path) {
  # create design
  temp_design <- data.frame(
    sampleID = sample_collection_metadata$sample_id,
    condition = sample_collection_metadata$tissue == tissue
  )

  temp_design["condition"][temp_design["condition"] == TRUE] <- tissue
  temp_design["condition"][temp_design["condition"] == FALSE] <- "other"
  # create switchlist
  temp_switchlist <- importRdata(
    isoformCountMatrix = merged_counts_iso,
    isoformRepExpression = cpm_iso,
    designMatrix = temp_design,
    isoformExonAnnoation = "/data/project/lasseigne_lab/TCH_scratch/data_ej/nextflow/bambu/extended_annotations.gtf",
    isoformNtFasta = "/data/project/lasseigne_lab/TCH_scratch/data_ej/gffread/isoform_sequences.fa",
    showProgress = FALSE
  )
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run satuRn
  switchlist_analyzed <- isoformSwitchTestSatuRn(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = TRUE
  )
  # rename object
  assign(paste0(tissue, "_switchlist_analyzed"),
    switchlist_analyzed,
    envir = .GlobalEnv
  )
  # save object
  saveRDS(switchlist_analyzed,
    paste0(save_path, "/", tissue, "_switchlist_saturn.Rds")
  )
}

# this function is for extracting significant genes from a switchlist
# save_path is typically here("results", "dtu_genes")
get_sig_genes <- function(tissue_switchlist, save_path) {
  # subset features
  sig_features <- dplyr::filter(
    tissue_switchlist$isoformFeatures,
    abs(dIF) > 0.1 & isoform_switch_q_value < 0.05
  )
  # get name
  name <- deparse(substitute(tissue_switchlist))
  name <- substr(name, 1, 4)
  # assign object
  assign(paste0(name, "_sig_features"),
    sig_features,
    envir = .GlobalEnv
  )
  # pull genes
  sig_genes <- unique(sig_features$gene_name)
  # assign object
  assign(paste0(name, "_sig_genes"),
    sig_genes,
    envir = .GlobalEnv
  )
  # save result
  write.table(sig_genes, 
              file = paste0(save_path, "/", name, "_others.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# this function is for plotting a single brain region/tissue object
# save_path is here("results", "plots", "satuRn_volcano")
create_volcano_plot_region <- function(tissue_obj, save_path) {
  # get name
  name <- deparse(substitute(tissue_obj))
  name <- substr(name, 1, 4)
  # plot
  volcano_plot <- ggplot(
    data = tissue_obj$isoformFeatures,
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
    ggtitle(paste0(name, " differentially used isoforms"))
  # save
  ggsave(
    paste0(
      save_path,
      "/", name, "_all_volcano.png"
    ),
    volcano_plot,
    width = 6, height = 4
  )
}

######################### 07_dtu_region_sex script #############################

# this function is for splitting out brain region counts with input "tissue"
split_region <- function(tissue) {
  # subset counts
  subset_counts <- subset(
    merged_counts_iso,
    select = c(
      "isoform_id",
      sample_collection_metadata$sample_id[
        sample_collection_metadata$tissue == tissue
      ]
    )
  )
  # subset cpm
  subset_cpm <- subset(
    cpm_iso,
    select = c(
      "isoform_id",
      sample_collection_metadata$sample_id[
        sample_collection_metadata$tissue == tissue
      ]
    )
  )
  # drop rows with 0 cpm
  subset_cpm <- subset_cpm[rowSums(subset_cpm[, -1]) != 0, ]
  # get name
  name <- substr(tissue, 1, 4)
  # assign objects
  assign(paste0(name, "_counts"), subset_counts, envir = .GlobalEnv)
  assign(paste0(name, "_cpm"), subset_cpm, envir = .GlobalEnv)
}

# this function is for making a switchlist object. it is sex-analysis specific.
make_switchlist <- function(tissue) {
  # get name
  name <- substr(tissue, 1, 4)
  # get objects for function
  assign("name_counts", get(paste0(name, "_counts")))
  assign("name_cpm", get(paste0(name, "_cpm")))
  # make design
  temp_design <- data.frame(
    sampleID = sample_collection_metadata$sample_id[
      sample_collection_metadata$tissue == tissue
    ],
    condition = sample_collection_metadata$sex[
      sample_collection_metadata$tissue == tissue
    ]
  )
  # make switchlist
  temp_switchlist <- importRdata(
    isoformCountMatrix = name_counts,
    isoformRepExpression = name_cpm,
    designMatrix = temp_design,
    isoformExonAnnoation = here(
      "data", "nextflow", "bambu", "extended_annotations.gtf"
    ),
    isoformNtFasta = here("data", "gffread", "isoform_sequences.fa"),
    showProgress = FALSE
  )
  # assign variable
  assign(paste0(name, "_sex_switchlist"), temp_switchlist, envir = .GlobalEnv)
}

# this function is for filtering data and runing satuRn for DTU analysis
# save path is typically here("data", "switchlist_objects") 
filter_run_saturn <- function(tissue, save_path) {
  # get name
  name <- substr(tissue, 1, 4)
  # assign name
  assign("temp_switchlist", get(paste0(name, "_sex_switchlist")))
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run satuRn
  switchlist_analyzed <- isoformSwitchTestSatuRn(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = TRUE
  )
if (!is.null(switchlist_analyzed)) {
  
  # examine summary
  extractSwitchSummary(switchlist_analyzed)
  
  # rename object
  assign(paste0(name, "_sex_switchlist_analyzed"),
         switchlist_analyzed,
         envir = .GlobalEnv
  )
  # save object
  saveRDS(switchlist_analyzed, paste0(save_path, "/",
                                      tissue, "_sex_switchlist_saturn.Rds")
  )
}
    
}

# this function is for sex split volcano plots for each brain region/tissue
# save_path is usually here("results", "plots", "satuRn_volcano")
create_sex_volcano_plot <- function(tissue_obj, save_path) {
  # get name
  name <- deparse(substitute(tissue_obj))
  name <- substr(name, 1, 4)
  # plot
  volcano_plot <- ggplot(
    data = tissue_obj$isoformFeatures,
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
      save_path,
      "/", name, "_sex_volcano.png"
    ),
    volcano_plot,
    width = 6, height = 4
  )
}

# this function is for getting significant isoforms from a certain region.
get_sig_isoforms <- function(region) {
  # get name of object
  name <- substr(region, 1, 4)
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

# this function is for getting significant genes and cutting off decimals
get_cut_sig_genes <- function(tissue) {
  # get name
  name <- substr(tissue, 1, 4)
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

# this function is for running gprofiler, skips plotting if no result returned
run_gprofiler <- function(tissue) {
  # get name
  name <- substr(tissue, 1, 4)
  # assign name
  assign(
    "temp_sex_sig_genes",
    get(paste0(name, "_sex_sig_genes"))
  )
  # run gprofiler
  temp_sex_gostres <- gost(
    query = temp_sex_sig_genes, organism = "mmusculus", ordered_query = FALSE,
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
    measure_underrepresentation = FALSE, evcodes = FALSE,
    user_threshold = 0.05, correction_method = "g_SCS", domain_scope = "custom",
    custom_bg = str_extract(
      unique(stri_sex_switchlist_analyzed$isoformFeatures$gene_id),
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

# this function is for getting gene symbols, though this function only works
# on a significant genes list, not the switchlist object.
# save path is here("results", "dtu_genes")
get_gene_symbols_sex <- function(tissue, save_path) {
  # get name of object
  name <- substr(tissue, 1, 4)
  # assign name
  assign("temp_sex_sig_genes", get(paste0(name, "_sex_sig_genes")))
  # get gene symbols from annotation dbi
  temp_sex_sig_gene_symbols <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = temp_sex_sig_genes,
    columns = c("SYMBOL", "GENENAME", "ENSEMBL"), keytype = "ENSEMBL"
  )
  # save results
  write.table(temp_sex_sig_gene_symbols$SYMBOL, 
              file = paste0(save_path, "/", name, "_sex.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  # rename object
  assign(paste0(name, "_sex_sig_gene_symbols"),
    temp_sex_sig_gene_symbols,
    envir = .GlobalEnv
  )
}

######################### 08_dtu_neuro_diseases script #########################

# this function is for converting human gene names to mouse gene names
convert_human_to_mouse <- function(gene_list){
  human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
                   host = "https://dec2021.archive.ensembl.org")
  mouse <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", 
                   host = "https://dec2021.archive.ensembl.org")
  # use biomaRt to get homologous genes
  genes <- getLDS(attributes = c("hgnc_symbol",'ensembl_gene_id'), 
                  filters = "hgnc_symbol", values = gene_list , mart = human, 
                  attributesL = c("mgi_symbol",'ensembl_gene_id'), 
                  martL = mouse, uniqueRows = TRUE)
  return(genes)
}

# this function is for pulling significant genes and getting overlap for a given
# brain region with three different gene lists.
compare_switching_genes <- function(brain_region) {
  # get name
  name <- brain_region
  # assign internal name
  assign("switchlist_analyzed", get(paste0(name, "_switchlist_analyzed")))
  # get significant genes
  sig_features <- dplyr::filter(
    switchlist_analyzed$isoformFeatures,
    abs(dIF) > 0.1 & isoform_switch_q_value < 0.05
  )
  # remove decimals from ENSEMBL ID
  sig_features$short_id <- str_extract(sig_features$gene_id,
                                       "ENSMUSG...........")
  # pull short gene ids
  dtu_genes <- unique(sig_features$short_id)
  
  # give name in global env
  assign(paste0(name, "_switching_genes"),
         dtu_genes,
         envir = .GlobalEnv
  )
  
  # compare to ad genes
  dtu_ad_genes <- intersect(dtu_genes, ad_mouse)
  
  # export genes
  assign(paste0(name, "_ad_genes"),
         dtu_ad_genes,
         envir = .GlobalEnv
  )
  
  # compare to psych genes
  dtu_psych_genes <- intersect(dtu_genes, psychiatric_mouse)
  # export genes
  assign(paste0(name, "_psych_genes"),
         dtu_psych_genes,
         envir = .GlobalEnv
  )
  
  # compare to CPAM genes
  dtu_cpam_genes <- intersect(dtu_genes, cpam_mouse)
  
  # export genes
  assign(paste0(name, "_cpam_genes"),
         dtu_cpam_genes,
         envir = .GlobalEnv
  )
  
}

######################### 09_dtu_isoform_switching script ######################

# this function is for adding and saving orfs for brain region
# save_path is here("data", "switchlist_objects", "orf_added")
add_save_orfs <- function(region, save_path) {
  # pull in switchlist
  assign("switchlist_analyzed", get(paste0(region, "_switchlist_analyzed")))
  # add open reading frames
  switchlist_analyzed <- addORFfromGTF(
    switchAnalyzeRlist = switchlist_analyzed,
    pathToGTF = here("data", "gencode_annotations",
                     "gencode.vM31.primary_assembly.annotation.gtf"))
  # add novel isoform orfs
  switchlist_analyzed <- analyzeNovelIsoformORF(
    switchlist_analyzed, analysisAllIsoformsWithoutORF = TRUE)
  # save
  saveRDS(switchlist_analyzed, paste0(save_path, "/",
      region,"_switchlist_orf.Rds"
    )
  )
  # assign object
  assign(paste0(region, "_switchlist_analyzed"), switchlist_analyzed,
         envir = .GlobalEnv)
}

# this function for adding and saving orfs for brain region (sex-specific)
# save_path is here("data", "switchlist_objects", "orf_added")
add_save_orfs_sex <- function(region, save_path) {
  # pull in switchlist
  assign("switchlist_analyzed", get(paste0(region, "_sex_switchlist_analyzed")))
  # add open reading frames
  switchlist_analyzed <- addORFfromGTF(
    switchAnalyzeRlist = switchlist_analyzed,
    pathToGTF = here("data", "gencode_annotations",
                     "gencode.vM31.primary_assembly.annotation.gtf"))
  # add novel isoform orfs
  switchlist_analyzed <- analyzeNovelIsoformORF(
    switchlist_analyzed, analysisAllIsoformsWithoutORF = TRUE)
  # save
  saveRDS(switchlist_analyzed, paste0(save_path, "/",
      region,"_sex_switchlist_orf.Rds"
    )
  )
  # assign object
  assign(paste0(region, "_sex_switchlist_analyzed"), switchlist_analyzed,
         envir = .GlobalEnv)
}