#' @import shiny
UploadSpectraUI <- function() {
  
  sidebarLayout(
    
    sidebarPanel = sidebarPanel(
      fluidRow(
        # Add radio buttons for user to specify the reference genome
        column(6, AddReferenceGenome2()),
        
        # Add radio buttons for user to specify the genomic region
        # from where the catalogs were generated
        column(6, AddRegionForSpectra())
      ),
      
      fluidRow(
        column(6, UploadSpectra()),
        
        column(6, 
               uiOutput(outputId = "showSpectraFromCatalog"), 
               br(),
               uiOutput(outputId = "sigAttributionFromCatalog"))
      ), # end fluidRown
      width = 4
    ), # end sidebarPanel
    
    mainPanel = mainPanel(
      fluidRow(
        column(6,
               wellPanel(
                 div(tags$b("Load example spectra", style = "color: #337ab7;")),
                 br(),
                 splitLayout(cellWidths = c("21%", "23%", "21%", "21%"),
                             actionButton(inputId = "preloadSBS96Spectra", 
                                          label = "SBS96"),
                             actionButton(inputId = "preloadSBS192Spectra", 
                                          label = "SBS192"),
                             actionButton(inputId = "preloadDBS78Spectra", 
                                          label = "DBS78"),
                             actionButton(inputId = "preloadIDSpectra", 
                                          label = "ID")),
                 br(),
                 downloadButton(outputId = "downloadSampleSpectra",
                                label = "Download example spectra"),
               ) # End wellPanel
        ) # end column
      )# end fluidRow
    ) # end mainPanel
  ) # End sidebarLayout
} # End function
