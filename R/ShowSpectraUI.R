#' @import shiny
ShowSpectraUI <- function() {
  fixedPage(
    shinyjs::useShinyjs(),
    sidebarLayout(
      
      sidebarPanel(
        uiOutput(outputId = "selectSampleFromUploadedVCF"),
        uiOutput(outputId = "selectSampleFromUploadedCatalog"),
        uiOutput(outputId = "buttonToSigAttribution")
      ),
      
      mainPanel(
        
        
        if (FALSE) {
          tabsetPanel(type = "tabs",
                      tabPanel("SBS96", plotOutput("SBS96plot")),
                      tabPanel("SBS192", plotOutput("SBS192plot")),
                      tabPanel("SBS1536", plotOutput("SBS1536plot")),
                      tabPanel("DBS78", plotOutput("DBS78plot")),
                      tabPanel("DBS136", plotOutput("DBS136plot")),
                      tabPanel("DBS144", plotOutput("DBS144plot")),
                      tabPanel("ID", plotOutput("IDplot"))
          )
        },
        
        uiOutput(outputId = "spectraPlotFromVCF"),
        uiOutput(outputId = "spectraPlotFromCatalog")
        
      )
      
    )
  )
}