library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)
library(ComplexHeatmap)
library(tidyverse)

# set ui
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  titlePanel("Visualizing Isoform Switches in Wild Type Mouse Brain"),
  # shiny app and paper overview tab
  tabsetPanel(
    id = "inTabset",
    tabPanel(
      "Welcome and About",
      h3("Welcome!"),
      p(
        "Welcome to the Shiny application created by Emma Jones in the",
        a("Lasseigne Lab", href = "https://www.lasseigne.org/"),
        "for visualizing long-read mouse brain RNA-sequencing data!"
      ),
      p(
        "To address differences in splicing across brain regions and sexes, we used long-read ONT RNA sequencing to sequence 40 mouse brain cDNA libraries from 10 mice and calculated differential transcript usage. We found that there is strong evidence of differential transcript usage across brain regions as well as differential expression at the gene and transcript level. We found that the brain region with the most differential expression and transcript usage is the cerebellum, potentially driven by differences in cell type composition. Overall, our findings suggest there is much differential splicing across brain regions and to a lesser extent, within brain regions across sexes. These differences in splicing could explain sex differences in prevalence and prognosis of various neurological and psychiatric disorders."
      ),
      img(
        src = "graphical_abstract.png",
        width = "825px",
        height = "225px"
      ),
      p(
        "If you are interested in checking out our work, feel free to check out our accompanying manuscript", a("here", href = "placeholder.com"),"."
      ),
      p(
        "This web application allows users to query their own genes of interest in our generated dataset and produce plots to explore their own hypotheses. Happy learning!"
        ),
      h3("What is differential transcript usage?"),
      p(
        "Differential transcript usage (DTU) measures changes in the proportions of transcripts expressed that a gene may have across groups."
      ),
      img(
        src = "DTU_example.png",
        width = "600px",
        height = "300px"),
      p(
        "In the figure above, there is an example on the left of differential gene expression across conditions, but no differential transcript usage since the transcript isoforms are still expressed in the same proportions. On the right, you can see that the transcript isoforms are expressed in different proportions, showing differential transcript usage without differential gene expression."
      ),
      p(
        "The isoform switch plots displayed are generated using the",
        a("IsoformSwitchAnalyzeR",
          href = "https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html"),
        "package created and maintained by Kristoffer Vitting-Seerup."
      ),
      img(src = "logo_only.png", width = "50px"),
      
      br(),
      br()
    ),
    
    # gene expression heatmap tab
    tabPanel(
      "Custom Gene Expression Heatmap",
      sidebarLayout(sidebarPanel(
        p(
          "Interested in gene-level expression (counts per million) of your genes of interest?"),
        p(
          "Use this tool to make a heatmap of you favorite genes (max recommended genes is around 30), and see how they cluster colored by tissue sample."
        ),
        selectizeInput(
          inputId = "gene_ids",
          label = "Input Gene Symbols or ENSEMBL IDs Here",
          choices = NULL,
          multiple = TRUE
        ),
        downloadButton("download_heatmap", "Download Raw Data"),
        downloadButton("download_heatmap_image", "Download Image", icon = icon("camera")),
        br(),
        br(),
        p(
          "Note: be mindful of expression values, heatmap color may be skewed when a very lowly expressed gene is plotted with more highly expressed genes."
        )
      ),
      mainPanel(imageOutput("heatmap_image")))
    ),
    
    # single brain region vs all others switchplot tab
    tabPanel(
      "Compare Single Brain Region to Others",
      sidebarLayout(
        sidebarPanel(
          p(
            "Are you interested in comparing differential transcript usage across a single brain region compared to other regions?"
          ),
          p(
            "If you are, select a gene name and a brain region to focus on. Alternatively, you can download the processed data to plot yourself!"
          ),
          selectizeInput(
            inputId = "gene_name1",
            label = "Input Gene Symbol Here",
            choices = NULL
          ),
          selectInput(
            inputId = "brain_region_1",
            label = "Select Brain Region",
            choices = list(
              "Cerebellum" = "cerebellum",
              "Cortex" = "cortex",
              "Hippocampus" = "hippocampus",
              "Striatum" = "striatum"
            ),
            selected = "cerebellum"
          ),
          downloadButton("download_data_1", "Download Raw Data"),
          downloadButton("download_image_1", "Download Image", icon = icon("camera")),
          br(),
          br(),
          p(
            "Note: if you do not see your gene of interest, it is likely it did not have at least 2 transcripts measured in our data or was not expressed in both samples and/or brain regions."
          )
        ),
        mainPanel(
          imageOutput("switchplot_1_image"))
      )
    ),
    
    # pairwise brain region switchplot tab
    tabPanel(
      "Compare Two Brain Regions",
      sidebarLayout(
        sidebarPanel(
          p(
            "Are you interested in comparing differential transcript usage pairwise across two specific brain regions?"
          ),
          p(
            "If you are, select a gene name and brain regions to focus on. Alternatively, you can download the processed data to plot yourself!"
          ),
          selectizeInput(
            inputId = "gene_name2",
            label = "Input Gene Symbol Here",
            choices = NULL
          ),
          selectInput(
            inputId = "condition1",
            label = "Select Brain Region 1",
            choices = list(
              "Cerebellum" = "cerebellum",
              "Cortex" = "cortex",
              "Hippocampus" = "hippocampus",
              "Striatum" = "striatum"
            ),
            selected = "cerebellum"
          ),
          selectInput(
            inputId = "condition2",
            label = "Select Brain Region 2",
            choices = list(
              "Cerebellum" = "cerebellum",
              "Cortex" = "cortex",
              "Hippocampus" = "hippocampus",
              "Striatum" = "striatum"
            ),
            selected = "cortex"
          ),
          downloadButton("download_data_2", "Download Raw Data"),
          downloadButton("download_image_2", "Download Image", icon = icon("camera")),
          br(),
          br(),
          p(
            "Note: if you do not see your gene of interest, it is likely it did not have at least 2 transcripts measured in our data or was not expressed in both samples and/or brain regions."
          )
        ),
        mainPanel(
          imageOutput("switchplot_2_image"))
      )
    ),
    
    # within brain region across sexes switchplot tab
    tabPanel(
      "Compare Sex Within a Brain Region",
      sidebarLayout(
        sidebarPanel(
          p(
            "Are you interested in comparing differential transcript usage within a single brain region across sexes?"
          ),
          p(
            "If you are, select a gene name and a brain region to focus on. Alternatively, you can download the processed data to plot yourself!"
          ),
          selectizeInput(
            inputId = "gene_name3",
            label = "Input Gene Symbol Here",
            choices = NULL
          ),
          selectInput(
            inputId = "brain_region_3",
            label = "Select Brain Region",
            choices = list(
              "Cerebellum" = "cerebellum",
              "Cortex" = "cortex",
              "Hippocampus" = "hippocampus",
              "Striatum" = "striatum"
            ),
            selected = "cerebellum"
          ),
          downloadButton("download_data_3", "Download Raw Data"),
          downloadButton("download_image_3", "Download Image", icon = icon("camera")),
          br(),
          br(),
          p(
            "Note: if you do not see your gene of interest, it is likely it did not have at least 2 transcripts measured in our data or was not expressed in both samples and/or brain regions."
          )
        ),
        
        mainPanel(imageOutput("switchplot_3_image"))
      )
    )
  )
)