#' @import shiny
xUploadSpectraUI <- function() {
  
  sidebarLayout(
    
    sidebarPanel = sidebarPanel(
      fluidRow(
        # Add radio buttons for user to specify the reference genome
        column(6, AddReferenceGenome2()),
        
        # Add radio buttons for user to specify the genomic region
        # from where the catalogs were generated
        column(6, AddRegion2()),
        
        column(6, UploadSpectra()),
        
        splitLayout(cellWidths = c("30%", "70%"),
                    uiOutput(outputId = "showSpectraFromCatalog"),
                    uiOutput(outputId = "sigAttributionsFromCatalog"))
        # column(6, offset = 6, uiOutput(outputId = "showSpectraFromCatalog"))
      )), # end fluidRow and sidebarPanel
    
    mainPanel = mainPanel(
      fluidRow(
        column(6,
               wellPanel(
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
                                          label = "ID")),
                 br(),
                 downloadButton(outputId = "downloadSampleSpectra",
                                label = "Download example spectra"),
               ) # End wellPanel
        ) # end column
      ) # end fluidRown
    ) # end mainPanel
  ) # End sidebarLayout
} # End function
