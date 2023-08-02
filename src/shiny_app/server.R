library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)

# global variables

region_all_switchlist_list_analyzed <- 
  readRDS("complete_switchlists/region_all_list_orf_de_pfam.Rds")

region_region_switchlist_analyzed <- 
  readRDS("complete_switchlists/region_region_orf_de_pfam.Rds")

region_sex_switchlist_list_analyzed <- 
  readRDS("complete_switchlists/region_sex_list_orf_de_pfam.Rds")

# code for server
function(input, output) {
  output$switchplot_1 <- renderPlot({
    
    switchPlot(region_all_switchlist_list_analyzed[[input$brain_region_1]],
               gene = input$gene_name1,
               plotTopology = FALSE)})
  
  output$switchplot_2 <- renderPlot({
    
    switchPlot(region_region_switchlist_analyzed,
               gene = input$gene_name2,
               condition1 = input$condition1,
               condition2 = input$condition2,
               plotTopology = FALSE)})
  
  output$switchplot_3 <- renderPlot({
    
    switchPlot(region_sex_switchlist_list_analyzed[[input$brain_region_3]],
               gene = input$gene_name3,
               plotTopology = FALSE)})
}