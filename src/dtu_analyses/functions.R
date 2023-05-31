# the purpose of this script is to include all functions used in dtu analyses
# in most functions, the input "region" is a character vector denoting brain region 

##################### pca_eda script ########################

# make plot function for PCA plot
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

##################### dtu_region_region script #####################

# this function is for filtering genes that are tissue x in condition 1 or 2
filter_genes_or <- function(tissue) {
  # subset genes of interest
  sig_genes_subset <- filter(
    sig_isoform_genes,
    condition_1 == tissue | condition_2 == tissue
  )
  sig_genes_subset <- unique(sig_genes_subset$gene_id)
  # name object
  name <- substr(tissue, 1, 4)
  assign(paste0("sig_isoform_genes_", name),
    sig_genes_subset,
    envir = .GlobalEnv
  )
}

# filtering genes for each comparison
# specifying tissue x for condition 1 and tissue y in condition 2
filter_genes_comp <- function(tissue1, tissue2) {
  # subset genes of interest
  sig_genes_subset <- filter(
    sig_isoform_genes,
    condition_1 == tissue1 & condition_2 == tissue2
  )
  sig_genes_subset <- unique(sig_genes_subset$gene_name)
  # name object
  name_1 <- substr(tissue1, 1, 4)
  name_2 <- substr(tissue2, 1, 4)
  assign(paste0("sig_isoform_genes_", name_1, "_", name_2),
         sig_genes_subset,
         envir = .GlobalEnv
  )
  write.table(sig_genes_subset, 
              file = here("results", "dtu_genes", 
                          paste0(name_1, "_", name_2, ".txt")
              ), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# create volcano plot function
# specifying tissue x for condition 1 and tissue y in condition 2
# save_path is here("results", "plots", "satuRn_volcano")
create_volcano_plot_comp <- function(tissue1, tissue2, save_path) {
  # create subset
  analyzed_subset <- dplyr::filter(
    region_region_switchlist_analyzed$isoformFeatures,
    condition_1 == tissue1 & condition_2 == tissue2
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
      x, " vs. ", y,
      " differentially used isoforms"
    ))
  # save
  ggsave(
    paste0(
      save_path,
      "/", x, "_", y, "_volcano.png"
    ),
    plot = volcano, width = 6, height = 4
  )
}

# create function
run_plot_gprofiler <- function(tissue_obj, save_path) {
  # get names
  name <- deparse(substitute(tissue_obj))
  name <- substr(name, nchar(name) - 9 + 1, nchar(name))
  # run gprofiler2
  gostres <- gost(
    query = tissue_obj,
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

##################### dtu_region_others script #####################

# create function for making switchlists and running saturn on tissue x
# nust specifiy save path as well
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
    isoformExonAnnoation = here(
      "data", "nextflow", "bambu", "extended_annotations.gtf"
    ),
    isoformNtFasta = here("data", "gffread", "isoform_sequences.fa"),
    showProgress = FALSE
  )
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run satuRn
  switchlist_analyzed <- isoformSwitchTestSatuRn(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = FALSE
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

# create function to get gene symbols for all genes of a tissue
# if you want this function to work for one of the "sex" conditions,
# you can just add the "_sex" as the input.
# save_path for this should be here("data", "switchlist_objects")
get_gene_symbols <- function(tissue, save_path) {
  # assign name
  assign("temp_switchlist", get(paste0(tissue, "_switchlist_analyzed")))
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
  assign(paste0(tissue, "_switchlist_analyzed"),
    temp_switchlist,
    envir = .GlobalEnv
  )
  saveRDS(temp_switchlist, 
    paste0(save_path, "/", tissue, "_switchlist_saturn.Rds")
  )
}
# make function for extracting significant genes
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

# create function for plotting a single brain region/tissue object
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

##################### dtu_region_sex script #####################

# make function to split out brain regions
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

# make function for making a switchlist
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

# make function to filter and run satuRn
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
# create function for sex split volcano plot for each brain region x
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
# make different function without reducing argument
# save_path is usually here("data", "switchlist_objects")
filter_run_saturn_noreduce <- function(tissue, save_path) {
  # get name
  name <- substr(tissue, 1, 4)
  # assign name
  assign("temp_switchlist", get(paste0(name, "_sex_switchlist")))
  # filter switchlist
  temp_switchlist <- preFilter(temp_switchlist, geneExpressionCutoff = NULL)
  # run saturn
  switchlist_analyzed <- isoformSwitchTestSatuRn(
    switchAnalyzeRlist = temp_switchlist,
    reduceToSwitchingGenes = FALSE
  )
  # rename object
  assign(paste0(name, "_sex_switchlist_analyzed"),
    switchlist_analyzed,
    envir = .GlobalEnv
  )
  # save object
  saveRDS(switchlist_analyzed, 
    paste0(save_path, "/", tissue, "_sex_switchlist_saturn.Rds")
  )
}

# make volcano plot function, but remove NAs
# save path is usually here("results", "plots", "satuRn_volcano")
create_sex_volcano_plot_rm <- function(tissue_obj, save_path) {
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

# write function to get significant isoforms
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

# make simple function to get significant genes and cutting off decimals
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
# make go analysis function that skips plotting if no result returned
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

# write function to get symbols for each list, though this function only works
# on the significant genes list, not the switchlist object.
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

###### dtu_neuro_diseases script #######
# write function to convert human to mouse gene names
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

# function for pulling significant genes and getting overlap for a given brain
# region with three different gene lists.
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

########## dtu_isoform_switching script #################
# function for adding and saving orfs for brain region x
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
  assign(paste0(x, "_switchlist_analyzed"), switchlist_analyzed,
         envir = .GlobalEnv)
}
# function for adding and saving orfs for brain region x (sex-specific)
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