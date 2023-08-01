library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)

# global variables

region_all_switchlist_list_analyzed <- 
  readRDS("complete_switchlists/region_all_list_orf_de_pfam.Rds")

#region_region_switchlist_analyzed <-

# code for server
function(input, output) {
  output$switchplot <- renderPlot({
    
    switchPlot(region_all_switchlist_list_analyzed[[input$brain_region]],
               gene = "Trio",
               plotTopology = FALSE)})
}