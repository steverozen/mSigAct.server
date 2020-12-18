#' @import shiny
SignatureAttributionUI <- function() {
  fixedPage(
    sidebarLayout(
      sidebarPanel(
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
        DT::dataTableOutput(outputId = "sigAetiologyTable")
      )
    )
  )
}