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
        "Welcome to the Shiny application created by Emma Jones for visualizing long-read mouse brain RNA-sequencing data!"
      ),
      p("The plots displayed are created from the",
        a("IsoformSwitchAnalyzeR", 
          href = "https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html"),
        "package."),
      actionButton(inputId = "start",
                   label = "Let's Get Started!")
    ),
    
    tabPanel(
      "Compare Single Brain Region",
      sidebarLayout(sidebarPanel(
        selectInput(
          inputId = "brain_region",
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
      
      mainPanel(plotOutput("switchplot")))
    ),
    tabPanel("Compare Double Brain Region"),
    tabPanel("Compare Single Region Across Sexes")
  )
)