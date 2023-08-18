library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)
library(ComplexHeatmap)
library(tidyverse)

# read in counts data
gene_exp_counts <- read.table("raw_counts/counts_gene.txt",
                              header = TRUE)

ensembl_ids <- str_extract(rownames(gene_exp_counts), "ENSMUSG...........")

# set ui
fluidPage(
  titlePanel("Visualizing Isoform Switches in Wild Type Mouse Brain"),
  
  # shiny app and paper overview tab
  tabsetPanel(
    id = "inTabset",
    tabPanel(
      "Welcome and About",
      br(),
      p(
        "Welcome to the Shiny application created by Emma Jones in the",
        a("Lasseigne Lab", href = "https://www.lasseigne.org/"),
        "for visualizing long-read mouse brain RNA-sequencing data!"
      ),
      img(src = "graphical_abstract.png",
          width = "825px",
          height = "225px"),
      p(
        "The plots displayed are created from the",
        a("IsoformSwitchAnalyzeR",
          href = "https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html"),
        "package created and maintained by Kristoffer Vitting-Seerup."
      ),
      actionButton(inputId = "start",
                   label = "Let's Get Started!"),
      img(src = "logo_only.png", width = "50px")
    ),
    
    # gene expression heatmap tab
    tabPanel(
      "Custom Gene Expression Heatmap",
      sidebarLayout(sidebarPanel(
        selectizeInput(
          inputId = "gene_ids",
          label = "Input ENSEMBL IDs Here",
          choices = NULL,
          multiple = TRUE
        )
      ),
      mainPanel(plotOutput("heatmap")))),
    
    # single brain region vs all others switchplot tab
    tabPanel(
      "Compare Single Brain Region",
      sidebarLayout(sidebarPanel(
        textInput(
          inputId = "gene_name1",
          label = "Gene Name",
          placeholder = "Enter Gene Name Here..."
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
        downloadButton("download_data_1", "Download Raw Data")
      ),
      
      mainPanel(plotOutput("switchplot_1")))
    ),
    
    # pairwise brain region switchplot tab
    tabPanel(
      "Compare Double Brain Region",
      sidebarLayout(sidebarPanel(
        textInput(
          inputId = "gene_name2",
          label = "Gene Name",
          placeholder = "Enter Gene Name Here..."
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
        downloadButton("download_data_2", "Download Raw Data")
      ),
      
      mainPanel(plotOutput("switchplot_2")))
    ),
    
    # within brain region across sexes switchplot tab
    tabPanel(
      "Compare Single Region Across Sexes",
      sidebarLayout(sidebarPanel(
        textInput(
          inputId = "gene_name3",
          label = "Gene Name",
          placeholder = "Enter Gene Name Here..."
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
        downloadButton("download_data_3", "Download Raw Data")
      ),
      
      mainPanel(plotOutput("switchplot_3")))
    )
  )
)