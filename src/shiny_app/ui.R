library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)
library(ComplexHeatmap)

fluidPage(
  titlePanel("Visualizing Isoform Switches in Wild Type Mouse Brain"),
  
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
    
    tabPanel(
      "Custom Gene Expression Heatmap",
      sidebarLayout(sidebarPanel(
        fileInput(
          inputId = "gene_ids",
          label = "Input Gene IDs Here",
          accept = c(".csv", ".tsv", ".txt")
        )
      ),
      mainPanel(plotOutput("heatmap")))),
    
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
        )
      ),
      
      mainPanel(plotOutput("switchplot_1")))
    ),
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
        )
      ),
      
      mainPanel(plotOutput("switchplot_2")))
    ),
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
        )
      ),
      
      mainPanel(plotOutput("switchplot_3")))
    )
  )
)