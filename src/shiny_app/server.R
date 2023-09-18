library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)
library(ComplexHeatmap)
library(tidyverse)
library(scales)
library(pals)

# global variables

# read in switchlists
region_all_switchlist_list_analyzed <-
  readRDS("complete_switchlists/region_all_list_orf_de_pfam.Rds")

region_region_switchlist_analyzed <-
  readRDS("complete_switchlists/region_region_orf_de_pfam.Rds")

region_sex_switchlist_list_analyzed <-
  readRDS("complete_switchlists/region_sex_list_orf_de_pfam.Rds")

# read in counts data
combined_cpm <- readRDS("data/combined_cpm.Rds")

combined_ids <- rownames(combined_cpm)

input <- list()

input$gene_ids <- head(rownames(combined_cpm))

# read in sample collection metadata
sample_collection_metadata <-
  readRDS("data/sample_collection_metadata.Rds")

# make tissue annotation colors constant instead of random
tissue_annotation <-
  HeatmapAnnotation(
    tissue = sample_collection_metadata$tissue,
    col = list(
      tissue = c(
        "cerebellum" = "#F8766D",
        "cortex" = "#7CAE00",
        "hippocampus" = "#00BFC4",
        "striatum" = "#C77CFF"
      )
    ),
    annotation_legend_param = list(
      tissue = list(
        nrow = 1
      )
    )
  )

# server module
server <- function(input, output, session) {
  # selectize input for heatmap
  updateSelectizeInput(session, "gene_ids",
    choices = combined_ids,
    server = TRUE
  )

  output$download_heatmap <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0("cpm_expression.csv")
    },
    content = function(file) {
      # write the data to the `file` that will be downloaded
      write.csv(combined_cpm, file)
    }
  )

  # download image
  output$download_heatmap_image <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0("custom_expression_heatmap.png")
    },
    content = function(file) {
      width <- session$clientData$output_heatmap_image_width
      height <- session$clientData$output_heatmap_image_height

      pixelratio <- session$clientData$pixelratio

      # generate the image file
      png(file,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )
      # plot code here
      plot <- Heatmap(
        as.matrix(combined_cpm[input$gene_ids, ]),
        name = "cpm",
        top_annotation = tissue_annotation,
        show_column_names = FALSE,
        col = ocean.deep(10),
        heatmap_legend_param = list(legend_direction = "horizontal")
      )

      # draw plot
      draw(plot,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom"
      )

      # turn off device
      dev.off()
    }
  )

  # create a dynamically-sized heatmap
  output$heatmap_image <- renderImage(
    {
      width <- session$clientData$output_heatmap_image_width
      height <- session$clientData$output_heatmap_image_height

      pixelratio <- session$clientData$pixelratio

      outfile <- tempfile(fileext = ".png")

      # generate the image file
      png(outfile,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )
      # plot code here
      plot <- Heatmap(
        as.matrix(combined_cpm[input$gene_ids, ]),
        name = "cpm",
        top_annotation = tissue_annotation,
        show_column_names = FALSE,
        col = ocean.deep(10),
        heatmap_legend_param = list(legend_direction = "horizontal")
      )

      # draw plot
      draw(plot,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom"
      )

      # turn off device
      dev.off()

      # return a list containing the filename
      list(
        src = outfile,
        width = width,
        height = height,
        alt = "A custom heatmap of gene expression values"
      )
    },
    deleteFile = TRUE
  )

  # selectize input for first switchplot tab
  updateSelectizeInput(
    session,
    "gene_name1",
    choices = region_region_switchlist_analyzed$isoformFeatures$gene_name,
    server = TRUE
  )

  # fetch data
  data_1 <-
    reactive(region_all_switchlist_list_analyzed[[
      input$brain_region_1
    ]]$isoformFeatures)

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

  # download image
  output$download_image_1 <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0(input$gene_name1, "_", input$brain_region_1, "_switchplot.png")
    },
    content = function(file) {
      width <- session$clientData$output_switchplot_1_image_width
      height <- session$clientData$output_switchplot_1_image_height

      pixelratio <- session$clientData$pixelratio

      # generate the image file
      png(file,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )

      # plot code here
      switchPlot(
        region_all_switchlist_list_analyzed[[input$brain_region_1]],
        gene = input$gene_name1,
        plotTopology = FALSE,
        localTheme = theme_bw(base_size = 11)
      )

      # turn off device
      dev.off()
    }
  )

  # a dynamically-sized switchplot for first tab
  output$switchplot_1_image <- renderImage(
    {
      width <- session$clientData$output_switchplot_1_image_width
      height <- session$clientData$output_switchplot_1_image_height

      pixelratio <- session$clientData$pixelratio

      outfile <- tempfile(fileext = ".png")

      # generate the image file
      png(outfile,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )

      # plot code here
      switchPlot(
        region_all_switchlist_list_analyzed[[input$brain_region_1]],
        gene = input$gene_name1,
        plotTopology = FALSE,
        localTheme = theme_bw(base_size = 11)
      )

      # turn off device
      dev.off()

      # return a list containing the filename
      list(
        src = outfile,
        width = width,
        height = height,
        alt = "A switch plot with gene of interest"
      )
    },
    deleteFile = TRUE
  )


  # selectize input for second switchplot tab
  updateSelectizeInput(
    session,
    "gene_name2",
    choices = region_region_switchlist_analyzed$isoformFeatures$gene_name,
    server = TRUE
  )

  # download data for second switchplot tab
  data_2 <- reactive({
    filter(
      region_region_switchlist_analyzed$isoformFeatures,
      condition_1 == input$condition1 &
        condition_2 == input$condition2
    )
  })

  output$download_data_2 <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0(
        input$condition1,
        "_",
        input$condition2,
        "_isoformFeatures.csv"
      )
    },
    content = function(file) {
      # write the data to the `file` that will be downloaded
      write.csv(data_2(), file)
    }
  )

  # download image
  output$download_image_2 <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0(
        input$gene_name2, "_", input$condition1, "_",
        input$condition2, "_switchplot.png"
      )
    },
    content = function(file) {
      width <- session$clientData$output_switchplot_2_image_width
      height <- session$clientData$output_switchplot_2_image_height

      pixelratio <- session$clientData$pixelratio

      # generate the image file
      png(file,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )

      # plot code here
      switchPlot(
        region_region_switchlist_analyzed,
        gene = input$gene_name2,
        condition1 = input$condition1,
        condition2 = input$condition2,
        plotTopology = FALSE,
        localTheme = theme_bw(base_size = 11)
      )

      # turn off device
      dev.off()
    }
  )

  # a dynamically-sized switchplot for second tab
  output$switchplot_2_image <- renderImage(
    {
      width <- session$clientData$output_switchplot_2_image_width
      height <- session$clientData$output_switchplot_2_image_height

      pixelratio <- session$clientData$pixelratio

      outfile <- tempfile(fileext = ".png")

      # generate the image file
      png(outfile,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )

      # plot code here
      switchPlot(
        region_region_switchlist_analyzed,
        gene = input$gene_name2,
        condition1 = input$condition1,
        condition2 = input$condition2,
        plotTopology = FALSE,
        localTheme = theme_bw(base_size = 11)
      )

      # turn off device
      dev.off()

      # return a list containing the filename
      list(
        src = outfile,
        width = width,
        height = height,
        alt = "A switch plot with gene of interest"
      )
    },
    deleteFile = TRUE
  )

  # selectize input for third switchplot tab
  updateSelectizeInput(
    session,
    "gene_name3",
    choices = region_region_switchlist_analyzed$isoformFeatures$gene_name,
    server = TRUE
  )
  # fetch data
  data_3 <-
    reactive(region_sex_switchlist_list_analyzed
    [[input$brain_region_3]]$isoformFeatures)

  # download data
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

  # download image
  output$download_image_3 <- downloadHandler(
    filename = function() {
      # set the suggested file name
      paste0(input$gene_name3, "_", input$brain_region_3, "_switchplot.png")
    },
    content = function(file) {
      width <- session$clientData$output_switchplot_3_image_width
      height <- session$clientData$output_switchplot_3_image_height

      pixelratio <- session$clientData$pixelratio

      # generate the image file
      png(file,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )

      # plot code here
      switchPlot(
        region_sex_switchlist_list_analyzed[[input$brain_region_3]],
        gene = input$gene_name3,
        plotTopology = FALSE,
        localTheme = theme_bw(base_size = 11)
      )

      # turn off device
      dev.off()
    }
  )

  # a dynamically-sized switchplot for third tab
  output$switchplot_3_image <- renderImage(
    {
      width <- session$clientData$output_switchplot_3_image_width
      height <- session$clientData$output_switchplot_3_image_height

      pixelratio <- session$clientData$pixelratio

      outfile <- tempfile(fileext = ".png")

      # generate the image file
      png(outfile,
        width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio
      )

      # plot code here
      switchPlot(
        region_sex_switchlist_list_analyzed[[input$brain_region_3]],
        gene = input$gene_name3,
        plotTopology = FALSE,
        localTheme = theme_bw(base_size = 11)
      )

      # turn off device
      dev.off()

      # return a list containing the filename
      list(
        src = outfile,
        width = width,
        height = height,
        alt = "A switch plot with gene of interest"
      )
    },
    deleteFile = TRUE
  )
}
