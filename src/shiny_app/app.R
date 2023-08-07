library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)
library(ComplexHeatmap)

source("ui.R")
source("server.R")

shinyApp(ui, server)