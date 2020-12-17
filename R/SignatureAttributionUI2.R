#' @import shiny
SignatureAttributionUI2 <- function() {
  fixedPage(
    sidebarLayout(
      sidebarPanel(
        #uiOutput(outputId = "selectSampleFromCatalogForAttribution2"),
        #uiOutput(outputId = "selectCancerType2"),
        
        #uiOutput(outputId = "chooseSigSubsetForSampleFromCatalog2"),
        #uiOutput(outputId = "addSig"),
        #br(),
        #uiOutput(outputId = "chooseMoreSigs"),
        #uiOutput(outputId = "analysisButton"),
        ######################################################
        uiOutput(outputId = "selectSampleForAttribution"),
        uiOutput(outputId = "selectCancerType"),
        
        uiOutput(outputId = "uploadedCatalogType"),
        uiOutput(outputId = "chooseCatalogType"),
        
        
        uiOutput(outputId = "chooseSigSubset"),
        uiOutput(outputId = "addSig"),
        br(),
        uiOutput(outputId = "chooseMoreSigs"),
        uiOutput(outputId = "analysisButton")
      ),
      
      mainPanel(
        # Must use DT::dataTableOutput instead of dataTableOutput,
        # otherwise, the data table will not show up
        #DT::dataTableOutput(outputId = "mytable"),
        DT::dataTableOutput(outputId = "mytable")
      )
    )
  )
}