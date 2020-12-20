#' @import shiny
ShowSpectraUI <- function() {
  fixedPage(
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