library(shiny)
library(IsoformSwitchAnalyzeR)
library(viridis)
library(ComplexHeatmap)
library(tidyverse)

source("ui.R")
source("server.R")

shinyApp(ui, server)