library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)

fluidPage(
  titlePanel("Visualizing Isoform Switches in Wild Type Mouse Brain"),
  
  tabsetPanel(
    tabPanel(
      "Welcome and About",
      br(),
      p(
        "Welcome to the Shiny application created by Emma Jones in the",
        a("Lasseigne Lab", href = "https://www.lasseigne.org/"),
        "for visualizing long-read mouse brain RNA-sequencing data!"
      ),
      p(
        "The plots displayed are created from the",
        a("IsoformSwitchAnalyzeR",
          href = "https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html"),
        "package created and maintained by Kristoffer Vitting-Seerup."
      ),
      actionButton(inputId = "start",
                   label = "Let's Get Started!")
    ),
    
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