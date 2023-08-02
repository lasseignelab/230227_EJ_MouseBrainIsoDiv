# Function script for use with the EJMouseBrainIsoDiv project and was written by
# Emma Jones. The purpose of this script is to include all functions used in DTU
# analyses. The input "tissue" or "region" is a character vector denoting brain
# region. This file is broken up into sections with functions for each
# individual numbered script or set of similar scripts.

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
        "%"
      ),
      y = paste0(
        deparse(substitute(secondPC)),
        ": ",
        round(var_explained[
          as.numeric(substr(deparse(substitute(secondPC)), 3, 4))
        ] * 100, 1),
        "%"
      )
    )

  if (length(unique(metadata[[color_len]])) > 10 &&
    color_class == "character" ||
    length(unique(metadata[[shape_len]])) > 10) {
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

######################### 05-07 dtu scripts ####################################
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
# save_path for this should be here("data", "switchlist_objects", "raw")
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

create_volcano_plot <-
  function(switchlist, condition1 = NULL, condition2 = NULL) {
    if (!is.null(condition1) && !is.null(condition2)) {
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

######################### 08_dtu_neuro_diseases script #########################

# this function is for converting human gene names to mouse gene names
convert_human_to_mouse <- function(gene_list) {
  human <- useMart("ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "https://dec2021.archive.ensembl.org"
  )
  mouse <- useMart("ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl",
    host = "https://dec2021.archive.ensembl.org"
  )
  # use biomaRt to get homologous genes
  genes <- getLDS(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters = "hgnc_symbol", values = gene_list, mart = human,
    attributesL = c("mgi_symbol", "ensembl_gene_id"),
    martL = mouse, uniqueRows = TRUE
  )
  return(genes)
}

# this function is for pulling significant genes and getting overlap for a given
# brain region with three different gene lists.
compare_switching_genes <- function(switchlist_index) {
  # pull chosen switchlist
  switchlist <- region_all_switchlist_analyzed[[switchlist_index]]
  #get name
  name <- names(region_all_switchlist_analyzed)[[switchlist_index]]
  # get significant genes
  sig_features <- dplyr::filter(
    switchlist$isoformFeatures,
    abs(dIF) > 0.1 & isoform_switch_q_value < 0.05
  )
  # remove decimals from ENSEMBL ID
  sig_features$short_id <- str_extract(
    sig_features$gene_id,
    "ENSMUSG..........."
  )
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

# this function is for adding and saving annotated and novel orfs
add_orfs <- function(switchlist_object) {
  # add open reading frames
  switchlist_analyzed <- addORFfromGTF(
    switchAnalyzeRlist = switchlist_object,
    pathToGTF = here(
      "data", "gencode_annotations",
      "gencode.vM31.primary_assembly.annotation.gtf"
    )
  )
  # add novel isoform orfs
  switchlist_analyzed <- analyzeNovelIsoformORF(
    switchlist_analyzed,
    analysisAllIsoformsWithoutORF = TRUE
  )
  # return object
  return(switchlist_analyzed)

}

# this function for adding and saving orfs for brain region (sex-specific)
# save_path is here("data", "switchlist_objects", "orf_added")
add_save_orfs_sex <- function(region, save_path) {
  # pull in switchlist
  assign("switchlist_analyzed", get(paste0(region, "_sex_switchlist_analyzed")))
  # add open reading frames
  switchlist_analyzed <- addORFfromGTF(
    switchAnalyzeRlist = switchlist_analyzed,
    pathToGTF = here(
      "data", "gencode_annotations",
      "gencode.vM31.primary_assembly.annotation.gtf"
    )
  )
  # add novel isoform orfs
  switchlist_analyzed <- analyzeNovelIsoformORF(
    switchlist_analyzed,
    analysisAllIsoformsWithoutORF = TRUE
  )
  # save
  saveRDS(switchlist_analyzed, paste0(
    save_path, "/",
    region, "_sex_switchlist_orf.Rds"
  ))
  # assign object
  assign(paste0(region, "_sex_switchlist_analyzed"), switchlist_analyzed,
    envir = .GlobalEnv
  )
}
