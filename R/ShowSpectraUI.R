#' @import shiny
ShowSpectraUI <- function() {
  fluidPage(
    shinyjs::useShinyjs(),
    sidebarLayout(
      
      sidebarPanel(
        uiOutput(outputId = "selectSampleFromUploadedVCF"),
        uiOutput(outputId = "selectSampleFromUploadedCatalog")
      ),
      
      mainPanel(
        uiOutput(outputId = "spectraPlotFromVCF"),
        uiOutput(outputId = "spectraPlotFromCatalog")
        
      )
      
    )
  )
}