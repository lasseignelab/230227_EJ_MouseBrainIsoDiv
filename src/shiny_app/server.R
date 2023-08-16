library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)
library(ComplexHeatmap)
library(tidyverse)

# global variables

# read in switchlists
region_all_switchlist_list_analyzed <- 
  readRDS("complete_switchlists/region_all_list_orf_de_pfam.Rds")

region_region_switchlist_analyzed <- 
  readRDS("complete_switchlists/region_region_orf_de_pfam.Rds")

region_sex_switchlist_list_analyzed <- 
  readRDS("complete_switchlists/region_sex_list_orf_de_pfam.Rds")

# read in counts data
gene_exp_counts <- read.table("raw_counts/counts_gene.txt",
                              header = TRUE)

ensembl_ids <- str_extract(rownames(gene_exp_counts), "ENSMUSG...........")

# read in sample collection metadata
sample_collection_metadata <- readRDS("raw_counts/sample_collection_metadata.RDS")

# calculate cpm
gene_exp_cpm <-
  do.call(cbind, lapply(seq_len(ncol(gene_exp_counts)), function(i) {
    gene_exp_counts[i] * 1e6 / sum(gene_exp_counts[i])
  }))

input <- list()

input$gene_ids <- head(rownames(gene_exp_cpm))

tissue_annotation <-
  HeatmapAnnotation(tissue = sample_collection_metadata$tissue)

# code for server
function(input, output, session) {
  
  observeEvent(input$start, {
    updateTabsetPanel(session, "inTabset", selected = "Custom Gene Expression Heatmap")
  })
  
  # make switchplot for first tab
  output$switchplot_1 <- renderPlot({
    
    switchPlot(region_all_switchlist_list_analyzed[[input$brain_region_1]],
               gene = input$gene_name1,
               plotTopology = FALSE)})
  
  # fetch data
  data_1 <- reactive({
    get(region_all_switchlist_list_analyzed[[input$brain_region_1]]$isoformFeatures)
  })
  
  output$download_data_1 <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0(input$brain_region_1, "_isoformFeatures.csv")
    },
    content = function(file) {
      # write the data to the `file` that will be downloaded
      write.csv(data_1(), file)
    }
  )
  
  # make switchplot for second tab
  output$switchplot_2 <- renderPlot({
    
    switchPlot(region_region_switchlist_analyzed,
               gene = input$gene_name2,
               condition1 = input$condition1,
               condition2 = input$condition2,
               plotTopology = FALSE)})
 
  # download data for second switchplot tab
  data_2 <- reactive({
    condition_subset <- filter(
      region_region_switchlist_analyzed$isoformFeatures,
      condition_1 == input$condition1 & condition_2 == input$condition2
    )
    get(condition_subset)
  })
  
  output$download_data_2 <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0(input$condition1, "_", input$condition2, "_isoformFeatures.csv")
    },
    content = function(file) {
      # write the data to the `file` that will be downloaded
      write.csv(data_2(), file)
    }
  )
  # make switchplot for third tab
  output$switchplot_3 <- renderPlot({
    
    switchPlot(region_sex_switchlist_list_analyzed[[input$brain_region_3]],
               gene = input$gene_name3,
               plotTopology = FALSE)})
  # fetch data
  data_3 <- reactive({
    get(region_sex_switchlist_list_analyzed[[input$brain_region_3]]$isoformFeatures)
  })
  
  output$download_data_3 <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0(input$brain_region_3, "_sex_isoformFeatures.csv")
    },
    content = function(file) {
      # write the data to the `file` that will be downloaded
      write.csv(data_3(), file)
    }
  )
  
  # selectize input for heatmap
  updateSelectizeInput(session, "gene_ids", choices = ensembl_ids, server = TRUE)
  
  # make heatmap
  output$heatmap <- renderPlot({
    
    Heatmap(as.matrix(gene_exp_cpm[input$gene_ids, ]), name = "cpm",
            top_annotation = tissue_annotation,
            show_column_names = FALSE)
    })
}