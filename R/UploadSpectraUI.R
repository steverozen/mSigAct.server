#' @import shiny
UploadSpectraUI <- function() {
  # List the first level UI elements here
  fixedPage(
    # Add the first row of control widgets
    fixedRow(
      # Add radio buttons for user to specify the reference genome
      column(6, AddReferenceGenome2()),
      
      # Add radio buttons for user to specify the genomic region
      # from where the catalogs were generated
      column(6, AddRegion2())
    ),
    
    # Add the next row of control widgets
    fixedRow(
      # Add a file upload control for user to upload spectra file
      column(6,
             div(tags$b("Load example spectra", style = "color: #337ab7;")),
             br(),
             splitLayout(cellWidths = c("20%", "20%", "20%", "20%"),
                         actionButton(inputId = "preloadSBS96Spectra", 
                                      label = "SBS96"),
                         actionButton(inputId = "preloadSBS192Spectra", 
                                      label = "SBS192"),
                         actionButton(inputId = "preloadDBS78Spectra", 
                                      label = "DBS78"),
                         actionButton(inputId = "preloadIDSpectra", 
                                      label = "ID"))
      )
    ),
    
    fixedRow(column(6, offset = 6,
                    UploadSpectra(),
                    downloadButton(outputId = "downloadSampleSpectra",
                                   label = "Download example spectra"),
    )),
    
    br(),
    
    # Add a download button for user to download sample spectra to test
    fixedRow(column(6, offset = 6,
                    splitLayout(cellWidths = c("30%", "70%"), 
                                uiOutput(outputId = "showSpectraFromCatalog"), 
                                uiOutput(outputId = "sigAttributionFromCatalog"))
    )),
    
    br(),
    
    fixedRow(column(6, offset = 6, 
                    uiOutput(outputId = "removeButton2")))
    
  )
}