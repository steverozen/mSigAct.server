#' @import shiny
SignatureAttributionUI2 <- function() {
  fixedPage(
    sidebarLayout(
      sidebarPanel(
        uiOutput(outputId = "selectSampleFromCatalogForAttribution2"),
        uiOutput(outputId = "selectCancerType2"),
        uiOutput(outputId = "uploadedCatalogType"),
        uiOutput(outputId = "chooseSigSubsetForSampleFromCatalog2"),
        uiOutput(outputId = "addSig"),
        br(),
        uiOutput(outputId = "chooseMoreSigs"),
        uiOutput(outputId = "analysisButton"),
        ######################################################
        uiOutput(outputId = "selectSampleFromVCFForAttribution"),
        uiOutput(outputId = "selectCancerTypeOfVCF"),
        uiOutput(outputId = "chooseCatalogType"),
        uiOutput(outputId = "chooseSigSubsetForVCF"),
        uiOutput(outputId = "addSigForVCF"),
        br(),
        uiOutput(outputId = "chooseMoreSigsForVCF"),
        uiOutput(outputId = "analysisButtonForVCF")
      ),
      
      mainPanel(
        # Must use DT::dataTableOutput instead of dataTableOutput,
        # otherwise, the data table will not show up
        DT::dataTableOutput(outputId = "mytable"),
        DT::dataTableOutput(outputId = "mytableForVCF")
      )
    )
  )
}