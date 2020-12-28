#' @import mSigAct
#' @import promises
#' @import ipc
#' @import shiny
#' @import shinydashboard
app_server <- function(input, output, session) {
  # List the first level callModules here
  tryCatch({
    # Predefine some values for later use
    fut <- NULL
    result_val <- reactiveVal()
    running <- reactiveVal(FALSE)
    interruptor <- ipc::AsyncInterruptor$new()
    future.planned <- FALSE
    
    # Create reactiveValues object
    # and set flag to 0 to prevent errors with adding NULL
    rv <- reactiveValues(catalogGeneratedFlag = 0)
    
    # Increase the file upload limit to 100MB
    options(shiny.maxRequestSize = 100*1024^2)
    
    # When user clicks the action link on Home page, direct user to the relevant tab
    observeEvent(input$linkToGenerateCatalogTab, {
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "generateCatalogTab")
    })
    
    observeEvent(input$linkToUploadSpectraTab, {
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "uploadSpectraTab")
    })
    
    # Create an empty list which can be used to store notification ids for
    # generating catalogs later
    ids <- list("error" = character(0), "warning" = character(0),
                "message" = character(0))
    
    # Create an empty list which can be used to store the return value of processing VCFs
    retval <- list()
    
    # Create variables which can be used to store values later
    catalog <- sig.universe <- zip.file.path <- NULL
    
    vcf.file.analyzed <- example.strelka.vcf.analyzed <- NULL
      
    example.mutect.vcf.analyzed <- NULL  
    
    list.of.catalogs <- catalog.path <- NA
    
    input.catalog.type <- reactiveVal(NULL)
    
    mutation.type <- NULL
    
    input.region <- reactiveVal(NULL)
      
    input.ref.genome <- reactiveVal(NULL)
    
    showSBS192Catalog <- TRUE
    
    attribution.results <- FALSE
    
    plotdata <- list(cossim = NULL, spect = NULL, 
                     reconstructed.catalog = NULL,
                     sig.universe = NULL, 
                     best.MAP.exp = NULL,
                     dat = NULL)
    
    dat <- data.frame(name = character(), spectrum = character(), 
                      proposed.aetiology = character())
    
    AdjustNumberOfCores <- 
      getFromNamespace(x = "Adj.mc.cores", ns = "mSigAct")
    
    #######################################################################
    # Reactive expressions
    #######################################################################
    choose.more.sigs <- reactive({
      more.sigs <- setdiff(colnames(total.signatures()), input$preselectedSigs)
      return(more.sigs)
    }) 
    
    total.signatures <- reactive(
      {
        if (is.null(input.catalog.type())) {
          return()
        }
        
        if (input.catalog.type() == "ID") {
          # The abudance information is not available for ID spectra, so 
          # we always use the human genome signatures for attribution analysis
          total.signatures <- COSMIC.v3.hg19.genome.ID.sigs
        } else {
          total.signatures <-
            COSMIC.v3.sigs[[input.ref.genome()]][[input.region()]][[input.catalog.type()]]
        }
        return(total.signatures)
      }
    )
    
    #######################################################################
    # Functions related to UploadVCFUI(), the first tab
    #######################################################################
    
    # Download sample VCFs when user clicks the button
    output$downloadsampleVCFs <- downloadHandler(
      filename = function() {
        "mSigAct-example-VCFs.zip"
      },
      content = function(file) {
        PrepareExampleVCFs(file)
      })
    
    
    observeEvent(input$runStrelkaVCFs, {
      # Create a temp path for the zip archive generated
      tmpdir <- tempfile()
      dir.create(tmpdir)
      zip.file.path <<- 
        file.path(tmpdir, "mSigAct-test-run-Strelka-VCFs-output.zip")
      
      results <- RunICAMSOnSampleStrelkaVCFs(session = session, output = output, 
                                             file = zip.file.path, ids = ids)
      list.of.catalogs <<- results$counts
      input.region("genome")
      input.ref.genome("hg19")
      
      # When the catalogs have been generated, increment rv$catalogGeneratedFlag
      rv$catalogGeneratedFlag <- rv$catalogGeneratedFlag + 1
      
      example.strelka.vcf.analyzed <<- TRUE
      example.mutect.vcf.analyzed <<- NULL
      vcf.file.analyzed <<- NULL
      
      shinyjs::hide(id = "clickToCreateCatalogs")
      shinyjs::hide(id = "spectraPlotFromVCF")
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
      
      output$downloadZipFile <- renderUI({
        downloadButton(outputId = "downloadCatalogs",
                       label = "Download catalogs")
      })
    })
    
    observeEvent(input$runMutectVCFs, {
      # Create a temp path for the zip archive generated
      tmpdir <- tempfile()
      dir.create(tmpdir)
      zip.file.path <<- 
        file.path(tmpdir, "mSigAct-test-run-Mutect-VCFs-output.zip")
      
      results <- RunICAMSOnSampleMutectVCFs(session = session, output = output, 
                                            file = zip.file.path, ids = ids)
      list.of.catalogs <<- results$counts
      input.region("genome")
      input.ref.genome("hg19")
      
      # When the catalogs have been generated, increment rv$catalogGeneratedFlag
      rv$catalogGeneratedFlag <- rv$catalogGeneratedFlag + 1
      
      example.mutect.vcf.analyzed <<- TRUE
      example.strelka.vcf.analyzed <<- NULL
      vcf.file.analyzed <<- NULL
      
      shinyjs::hide(id = "clickToCreateCatalogs")  
      shinyjs::hide(id = "spectraPlotFromVCF")
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
      
      output$downloadZipFile <- renderUI({
        downloadButton(outputId = "downloadCatalogs",
                       label = "Download catalogs")
      })
    })
    
    # When user uploads VCF files, then show the action button "Create catalogs"
    observeEvent(input$vcf.files, {
      output$clickToCreateCatalogs  <- renderUI({
        actionButton(inputId = "createCatalogs",
                     label = "Create catalogs",
                     style = "color: #fff; background-color: #337ab7;
                             border-color: #2e6da4;")
      })
      
      shinyjs::show(id = "clickToCreateCatalogs")  
      shinyjs::hide(id = "spectraPlotFromVCF")
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
    }) 
    
    # When user clicks "Create catalogs" button, generate catalogs and show
    # "Download catalogs" button
    observeEvent(input$createCatalogs, {
      
      # Check the arguments for generating catalogs from VCF
      errors <- CheckInputsForVCF(input)
      ids$error <<- append(ids$error, AddErrorMessage(errors))
      
      req(input$variantCaller, input$ref.genome, input$region, input$vcf.files)
      if (input$variantCaller == "unknown") {
        if (is.null(input$mergeSBS)) {
          return()
        }
      }
      input.ref.genome(input$ref.genome)
      input.region(input$region)
      
      # Create a temp path for the zip archive generated
      tmpdir <- tempfile()
      dir.create(tmpdir)
      zip.file.path <<- file.path(tmpdir, "mSigAct-catalogs-output.zip")
      
      result <- ProcessVCFs(input = input, output = output, 
                            file = zip.file.path, ids = ids)
      retval <<- result$retval
      old.error.ids <- ids$error
      ids <<- result$ids
      ids$error <- append(ids$error, old.error.ids)
      
      list.of.catalogs <<- retval$counts
      #density.catalog <- retval$density
      
      # When catalogs have been generated, increment rv$catalogGeneratedFlag
      rv$catalogGeneratedFlag <- rv$catalogGeneratedFlag + 1
      
      vcf.file.analyzed <<- TRUE
      example.strelka.vcf.analyzed <<- NULL
      example.mutect.vcf.analyzed <<- NULL
      
      output$downloadZipFile <- renderUI({
        downloadButton(outputId = "downloadCatalogs",
                       label = "Download catalogs")
      })
      
      shinyjs::show(id = "downloadZipFile")
      shinyjs::show(id = "selectSampleFromUploadedVCF")
    })
    
    # When user clicks "Download catalogs" button, prepare the zip file for
    # downloading
    output$downloadCatalogs <- downloadHandler(
      filename = function() {
        if (is.null(vcf.file.analyzed)) {
          if (is.null(example.strelka.vcf.analyzed)) {
            return("mSigAct-test-run-Mutect-VCFs-output.zip")
          } else {
            return("mSigAct-test-run-Strelka-VCFs-output.zip")
          }
        } else {
          return("mSigAct-catalogs-output.zip")
        } 
      },
      content = function(file) {
        file.copy(from = zip.file.path, to = file)
      })
    
    # When user uploads new VCF, hide the "Download catalogs",  "Show spectra"
    # and "Signature attribution" button
    observeEvent(input$vcf.files, {
      output$showSpectraFromVCF <- NULL
      output$sigAttributionFromVCF <- NULL
      output$downloadZipFile <- NULL
      shinyjs::hide(id = "showSpectraFromVCF")
      shinyjs::hide(id = "sigAttributionFromVCF")
      shinyjs::hide(id = "downloadZipFile")
    })
    
    # When user has generated catalog from VCF, then show the "Show spectra" and
    # "Signature attribution" button
    observeEvent(rv$catalogGeneratedFlag, {
      # Make the input.catalog.type value to be NULL
      # After newly generating catalogs from VCF, wait for the user to select
      # catalog type for attribution analysis
      input.catalog.type(NULL)
      
      output$showSpectraFromVCF <- renderUI({
        actionButton(inputId = "showSpectraOfVCF", label = "Show spectra",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      })
      
      output$sigAttributionFromVCF <- renderUI({
        actionButton(inputId = "sigAttributionOfVCF", label = "Signature attribution",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      })
      
      shinyjs::show(id = "showSpectraFromVCF")
      shinyjs::show(id = "sigAttributionFromVCF")
      
      output$selectSampleFromUploadedVCF <- renderUI(
        {
          sample.names <- colnames(list.of.catalogs[[1]])
          
          tagList(
            radioButtons(inputId = "sampleNameFromUploadedVCF",
                         label = "Select sample",
                         choices = sample.names,
                         selected = character(0)),
            
            actionButton(inputId = "clickToSigAttributionForVCF",
                         label = "Signature attribution",
                         style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
          )
        }
      )
      #shinyjs::show(id = "sampleNameFromUploadedVCF")
    }, ignoreInit = TRUE)
    
    observeEvent(input$showSpectraOfVCF, {
      # Must delete the previous plot and widget before hiding them
      output$spectraPlotFromCatalog <- NULL
      output$selectSampleFromUploadedCatalog <- NULL
      
      # Hide the previous plot from uploaded catalog
      shinyjs::hide(id = "spectraPlotFromCatalog")
      shinyjs::hide(id = "selectSampleFromUploadedCatalog ")
      
      shinyjs::show(id = "selectSampleFromUploadedVCF")
      
      shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "showSpectraTab")
    })
    
    observeEvent(input$sigAttributionOfVCF, {
      # Must delete the previous plot and widget before hiding them
      output$spectraPlotFromCatalog <- NULL
      output$selectSampleFromUploadedCatalog <- NULL
      
      # Hide the previous plot from uploaded catalog
      shinyjs::hide(id = "spectraPlotFromCatalog")
      shinyjs::hide(id = "selectSampleFromUploadedCatalog ")
      
      shinyjs::show(id = "selectSampleFromUploadedVCF")
      
      shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab2")
    })
    
    #######################################################################
    # Start of functions related to ShowSpectraUI
    #######################################################################
    
    # Download sample catalogs when user clicks the button
    output$downloadExampleSpectra <- downloadHandler(
      filename = function() {
        "mSigAct-example-spectra.zip"
      },
      content = function(file) {
        PrepareExampleSpectra(file)
      })
    
    # When user selects the sample from uploaded VCF, show
    # the sample's mutational spectrum
    observeEvent(input$sampleNameFromUploadedVCF, {
      output$spectraPlotFromVCF <- renderUI (
        {
          output$SBS96plot <- renderPlot({
            catSBS96 <-
              list.of.catalogs$catSBS96[, input$sampleNameFromUploadedVCF, drop = FALSE]
            ICAMS::PlotCatalog(catSBS96)
          }, height = 230, width = 800)
          
          
          output$SBS192plot <- renderPlot({
            catSBS192 <-
              list.of.catalogs$catSBS192[, input$sampleNameFromUploadedVCF, drop = FALSE]
            ICAMS::PlotCatalog(catSBS192)
          }, height = 250, width = 800)
          
          if (FALSE) {
            output$SBS12plot <- renderPlot({
              catSBS192 <-
                list.of.catalogs$catSBS192[, input$sampleNameFromUploadedVCF, drop = FALSE]
              ICAMS::PlotCatalog(catSBS192, plot.SBS12 = TRUE)
            }, height = 350, width = 350)
            
          }
           
          output$SBS1536plot <- renderPlot({
            catSBS1536 <-
              list.of.catalogs$catSBS1536[, input$sampleNameFromUploadedVCF, drop = FALSE]
            ICAMS::PlotCatalog(catSBS1536)
          }, height = 800, width = 800)
          
          output$DBS78plot <- renderPlot({
            catDBS78 <-
              list.of.catalogs$catDBS78[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catDBS78)
          }, height = 250, width = 800)
          
          output$DBS136plot <- renderPlot({
            catDBS136 <-
              list.of.catalogs$catDBS136[, input$sampleNameFromUploadedVCF, drop = FALSE]
            ICAMS::PlotCatalog(catDBS136)
          }, height = 500, width = 700)
          
          output$DBS144plot <- renderPlot({
            catDBS144 <-
              list.of.catalogs$catDBS144[, input$sampleNameFromUploadedVCF, drop = FALSE]
            ICAMS::PlotCatalog(catDBS144)
          }, height = 350, width = 350)
          
          output$IDplot <- renderPlot({
            catID <-
              list.of.catalogs$catID[, input$sampleNameFromUploadedVCF, drop = FALSE]
            ICAMS::PlotCatalog(catID)
          }, height = 230, width = 800)
          
          tabsetPanel(type = "tabs",
                      tabPanel("SBS96", plotOutput("SBS96plot")),
                      tabPanel("SBS192", plotOutput("SBS192plot")),
                      #tabPanel("SBS12", plotOutput("SBS12plot")),
                      tabPanel("SBS1536", plotOutput("SBS1536plot")),
                      tabPanel("DBS78", plotOutput("DBS78plot")),
                      tabPanel("DBS136", plotOutput("DBS136plot")),
                      tabPanel("DBS144", plotOutput("DBS144plot")),
                      tabPanel("ID", plotOutput("IDplot"))
          )
        })
      shinyjs::show(id = "spectraPlotFromVCF")
    })
    
    # When user selects the sample from uploaded catalog, show
    # the sample's mutational spectrum
    observeEvent(input$selectedSampleFromUploadedCatalog, {
      output$SBS96plot <- renderPlot({
        catSBS96 <-
          list.of.catalogs$catSBS96[, input$sampleNameFromUploadedVCF, drop = FALSE]
        ICAMS::PlotCatalog(catSBS96)
      }, height = 230, width = 800)
      
      
      output$SBS192plot <- renderPlot({
        catSBS192 <-
          list.of.catalogs$catSBS192[, input$sampleNameFromUploadedVCF, drop = FALSE]
        ICAMS::PlotCatalog(catSBS192)
      }, height = 250, width = 800)
      
      output$SBS12plot <- renderPlot({
        catSBS192 <-
          list.of.catalogs$catSBS192[, input$sampleNameFromUploadedVCF, drop = FALSE]
        ICAMS::PlotCatalog(catSBS192, plot.SBS12 = TRUE)
      }, height = 350, width = 350)
      
      output$SBS1536plot <- renderPlot({
        catSBS1536 <-
          list.of.catalogs$catSBS1536[, input$sampleNameFromUploadedVCF, drop = FALSE]
        ICAMS::PlotCatalog(catSBS1536)
      }, height = 800, width = 800)
      
      output$DBS78plot <- renderPlot({
        catDBS78 <-
          list.of.catalogs$catDBS78[, input$sampleNameFromUploadedVCF, drop = FALSE]
        PlotCatalog(catDBS78)
      }, height = 250)
      
      output$DBS136plot <- renderPlot({
        catDBS136 <-
          list.of.catalogs$catDBS136[, input$sampleNameFromUploadedVCF, drop = FALSE]
        ICAMS::PlotCatalog(catDBS136)
      }, height = 500, width = 700)
      
      output$DBS144plot <- renderPlot({
        catDBS144 <-
          list.of.catalogs$catDBS144[, input$sampleNameFromUploadedVCF, drop = FALSE]
        ICAMS::PlotCatalog(catDBS144)
      }, height = 350, width = 350)
      
      output$IDplot <- renderPlot({
        catID <-
          list.of.catalogs$catID[, input$sampleNameFromUploadedVCF, drop = FALSE]
        ICAMS::PlotCatalog(catID)
      }, height = 230, width = 800)
      
      
      if (input.catalog.type() == "SBS96") {
        height <- 230
        width <- 800
      } else if (input.catalog.type() == "SBS192") {
        height <- 250
        width <- 800
      } else if (input.catalog.type() == "SBS1536") {
        height <- 800
        width <- 800
      } else if (input.catalog.type() == "DBS78") {
        height <- 250
        width <- 800
      } else if (input.catalog.type() == "DBS136") {
        height <- 500
        width <- 700
      } else if (input.catalog.type() == "DBS144") {
        height <- 350
        width <- 350
      } else if (input.catalog.type() == "ID") {
        height <- 230
        width <- 800
      } 
      
      output$spectraPlotFromCatalog <- renderUI(
        {
          output$spectrum <- renderPlot(
            {
              PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                                  drop = FALSE])
            }, width = width, height = height)
          plotOutput(outputId = "spectrum")
        }
      )
      shinyjs::show(id = "spectraPlotFromCatalog")
    })
    
    # After showing the spectra plot from catalog, show an actionButton to redirect
    # to "Signature attribution" tab if the spectra catalog type is supported for
    # signature attribution
    observeEvent(CheckArgumentsForSpectra(), {
      
      req(input$ref.genome2, input$region2)
      input.ref.genome(input$ref.genome2)
      input.region(input$region2)
      
      # Delete the previous spectra plot
      output$spectraPlotFromCatalog <- NULL
      
      if (showSBS192Catalog == FALSE) {
        return()
      }
      
      output$selectSampleFromUploadedCatalog <-
        renderUI(
          {
            # Change the region information for SBS192 catalog
            my.region <- ChangeRegionForSBS192Catalog(input, catalog.path)
            
            catalog <<- ICAMS::ReadCatalog(file = catalog.path,
                                           ref.genome = input$ref.genome2,
                                           region = my.region)
            input.catalog.type(CheckCatalogType(catalog)) 
            mutation.type <<- CheckCatalogType(catalog)
            
            sample.names <- colnames(catalog)
            
            tagList(
              radioButtons(inputId = "selectedSampleFromUploadedCatalog",
                           label = "Select spectrum to view",
                           choices = sample.names,
                           selected = character(0)),
              
              if (input.catalog.type() %in% c("SBS96", "SBS192", "DBS78", "ID")) {
                actionButton(inputId = "clickToSigAttributionForSpectra",
                             label = "Signature attribution",
                             style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
              }
              
            )
          }
        ) # end of renderUI
    })
    
    # When user clicks the action button on Show spectra page, direct user to the relevant tab
    observeEvent(input$clickToSigAttributionForSpectra, {
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab2")
    })
    
    observeEvent(input$clickToSigAttributionForVCF, {
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab2")
    })
    
    # When user submit new catalog for analysis, remove the previous plots
    observeEvent(input$showSpectraOfCatalog, {
      output$SBS96plot <- NULL
      output$SBS192plot <- NULL
      output$SBS1536plot <- NULL
      output$DBS78plot <- NULL
      output$DBS136plot <- NULL
      output$DBS144plot <- NULL
      output$IDplot <- NULL
    })
    
    ########################################################################
    # Start of functions related to UploadSpectraUI
    ########################################################################
    ShowTwoButtons <- function() {
      output$showSpectraFromCatalog <- renderUI(
        actionButton(inputId = "showSpectraOfCatalog", label = "Show spectra",
                     style="color: #fff; background-color: #337ab7;
                          border-color: #2e6da4;")
      )
      output$sigAttributionFromCatalog <- renderUI(
        actionButton(inputId = "sigAttributionOfCatalog", 
                     label = "Signature attribution",
                     style="color: #fff; background-color: #337ab7;
                          border-color: #2e6da4;"))
    }
    
    observeEvent(input$preloadSBS96Spectra, {
      
      # Hide the UI element previously triggered by uploaded VCF
      shinyjs::hide(id = "chooseCatalogType")
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "genome")
      catalog.path <<- system.file("extdata/SBS96-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type("SBS96") 
      mutation.type <<- "SBS96"
      ShowTwoButtons()
      
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
    })
    
    observeEvent(input$preloadSBS192Spectra, {
      # Hide the UI element previously triggered by uploaded VCF
      shinyjs::hide(id = "chooseCatalogType")
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "genome")
      catalog.path <<- system.file("extdata/SBS192-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type("SBS192") 
      mutation.type <<- "SBS192"
      ShowTwoButtons()
      
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
    })
    
    observeEvent(input$preloadDBS78Spectra, {
      # Hide the UI element previously triggered by uploaded VCF
      shinyjs::hide(id = "chooseCatalogType")
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "genome")
      catalog.path <<- system.file("extdata/DBS78-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type("DBS78") 
      mutation.type <<- "DBS78"
      ShowTwoButtons()
      
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
    })
    
    observeEvent(input$preloadIDSpectra, {
      # Hide the UI element previously triggered by uploaded VCF
      shinyjs::hide(id = "chooseCatalogType")
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "genome")
      catalog.path <<- system.file("extdata/ID-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type("ID")
      mutation.type <<- "ID"
      ShowTwoButtons()
      
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
    })
    
    observeEvent(input$upload.spectra, {
      # catalog.info is a data frame that contains one row for each uploaded file,
      # and four columns "name", "size", "type" and "datapath".
      # "name": The filename provided by the web browser.
      # "size": The size of the uploaded data, in bytes.
      # "type": The MIME type reported by the browser.
      # "datapath": The path to a temp file that contains the data that was uploaded.
      catalog.info <- input$upload.spectra
      catalog.path <<- catalog.info$datapath
      ShowTwoButtons()
      
      shinyjs::hide(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::hide(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinyjs::hide(selector = '#panels li a[data-value=attributionResultsTab]')
    })
    
    # Check the arguments for uploaded spectra
    observeEvent(
      CheckArgumentsForSpectra(),
      {
        if (TwoActionButtonsClicked(input) == FALSE) {
          return()
        } else {
          retval <- CheckInputsForSpectra(input, catalog.path)
          ids$error <<- append(ids$error, AddErrorMessage(retval$error))
          showSBS192Catalog <<- retval$SBS192.check
        }
      }
    )
    
    # When user clicks the action button on "Upload Spectra" page, direct user to
    # the relevant tab
    observeEvent(input$showSpectraOfCatalog, {
      if (TwoActionButtonsClicked(input) == FALSE) {
        return()
      } else if (showSBS192Catalog == FALSE) {
        return()
      } else {
        req(input$ref.genome2, input$region2)
        input.ref.genome(input$ref.genome2)
        input.region(input$region2)
        
        # Change the region information for SBS192 catalog
        my.region <- ChangeRegionForSBS192Catalog(input, catalog.path)
        
        catalog <<- 
          ICAMS::ReadCatalog(file = catalog.path, ref.genome = input$ref.genome2,
                             region = my.region)
        input.catalog.type(CheckCatalogType(catalog))
        mutation.type <<- CheckCatalogType(catalog)
        
        # Hide the widgets for previous uploaded VCF
        output$spectraPlotFromVCF <- NULL
        output$selectSampleFromUploadedVCF <- NULL
        shinyjs::hide(id = "spectraPlotFromVCF")
        shinyjs::hide(id = "selectSampleFromUploadedVCF")
        
        shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
        
        if (input.catalog.type() %in% c("SBS96", "SBS192", "DBS78", "ID")) {
          shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
        }
        shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                       selected = "showSpectraTab")
      }
    })
    
    observeEvent(input$sigAttributionOfCatalog, {
      
      if (TwoActionButtonsClicked(input) == FALSE) {
        return()
      } else if (showSBS192Catalog == FALSE) {
        return()
      } else {
        req(input$ref.genome2, input$region2)
        input.ref.genome(input$ref.genome2)
        input.region(input$region2)
        
        # Change the region information for SBS192 catalog
        my.region <- ChangeRegionForSBS192Catalog(input, catalog.path)
        
        catalog <<- 
          ICAMS::ReadCatalog(file = catalog.path, ref.genome = input$ref.genome2,
                             region = my.region)
        input.catalog.type(CheckCatalogType(catalog))
        mutation.type <<- CheckCatalogType(catalog)
        
        if (!input.catalog.type() %in% c("SBS96", "SBS192", "DBS78", "ID")) {
          showNotification(ui = "Error:", 
                           action = paste0("Can only do signature attribution ", 
                                           "for SBS96, SBS192, DBS78 and ID"),
                           type = "error", duration = NULL)
          return()
        }
        # Hide the widgets for previous uploaded VCF
        output$spectraPlotFromVCF <- NULL
        output$selectSampleFromUploadedVCF <- NULL
        shinyjs::hide(id = "spectraPlotFromVCF")
        shinyjs::hide(id = "selectSampleFromUploadedVCF")
        
        
        output$uploadedCatalogType <- renderUI(
          {
            p(tags$b(paste0("Mutation type: ", mutation.type)))
          }
        )
        shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
        shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
        shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                       selected = "sigAttributionTab2")
      }
    })
    
    ########################################################################
    # Start of functions related to SignatureAttributionUI
    ########################################################################
    # When user clicks either of the two actionButtons "Show spectra", 
    # "Signature attribution" on UploadSpectraUI page, then update
    # the shiny widgets used for attribution analysis
    CheckArgumentsForSpectra <- reactive({
      list(input$showSpectraOfCatalog, input$sigAttributionOfCatalog)
    })
    
    TwoActionButtonsClicked <- function(input) {
      if(is.null(input$showSpectraOfCatalog) && is.null(input$sigAttributionOfCatalog)){
        return(FALSE)
      }
      if(input$showSpectraOfCatalog == 0 && input$sigAttributionOfCatalog == 0){
        return(FALSE)
      }
      return(TRUE)
    }
    
    observeEvent(
      CheckArgumentsForSpectra(),
      { 
        
        req(input$ref.genome2, input$region2)
        input.ref.genome(input$ref.genome2)
        input.region(input$region2)
        
        # Check whether the SBS192 catalog has the correct "region" parameter
        if (showSBS192Catalog == FALSE) {
          return()
        }
        
        output$selectSampleForAttribution <- renderUI(
          { 
            # Change the region information for SBS192 catalog
            my.region <- ChangeRegionForSBS192Catalog(input, catalog.path)
            catalog <<- ICAMS::ReadCatalog(file = catalog.path,
                                           ref.genome = input$ref.genome2,
                                           region = my.region)
            input.catalog.type(CheckCatalogType(catalog))
            mutation.type <<- CheckCatalogType(catalog)
            sample.names <- colnames(catalog)
            selectInput(inputId = "selectedSampleForAttribution",
                        label = "Select sample",
                        choices = sample.names)
          }
        )
        
        output$selectCancerType <- renderUI(
          {
            cancer.types <-
              c("Unknown", colnames(CancerTypeToExposureStatData()))
            selectInput(inputId = "selectedCancerType",
                        label = "Select cancer type",
                        choices = cancer.types)
          }
        )
        
        # Hide the UI element previously generated by uploading VCFs
        output$chooseCatalogType <- NULL
        shinyjs::hide(id = "chooseCatalogType")
        
        
        # Hide the previous signature aetiology table
        output$sigAetiologyTable <- NULL
        shinyjs::hide(id = "sigAetiologyTable")
        
        
        output$uploadedCatalogType <- renderUI(
          {
            p(tags$b(paste0("Mutation type: ", mutation.type)))
          }
        )
        
        shinyjs::show(id = "uploadedCatalogType")
      })
    
    observeEvent(input$selectedCancerType, {
      if (is.null(input.catalog.type())) {
        return()
      }
      
      if (!input.catalog.type() %in% c("SBS96", "SBS192", "DBS78", "ID")) {
        showNotification(ui = "Error:", 
                         action = paste0("Can only do signature attribution ", 
                                         "for SBS96, SBS192, DBS78 and ID"),
                         type = "error", duration = NULL)
        return()
      }
      
      if (input$selectedCancerType != "Unknown") {
        ShowPreselectedSigs(input, output,input.catalog.type())
        } # end of if statement
    
      # Show the actionButton for user to add more signatures
      output$addSig <- renderUI(
        actionButton(inputId = "addMoreSigs", label = "Add more signatures",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      )
    }, ignoreInit = TRUE) # end of observeEvent
    
    
    # When user clicks either of the two actionButtons "Show spectra", 
    # "Signature attribution" on UploadVCFUI page, then update
    # the shiny widgets used for attribution analysis
    CheckArgumentsForVCF <- reactive({
      list(input$showSpectraOfVCF, input$sigAttributionOfVCF)
    })
    
    TwoVCFActionButtonsClicked <- function(input) {
      if(is.null(input$showSpectraOfVCF) && is.null(input$sigAttributionOfVCF)){
        return(FALSE)
      }
      if(input$showSpectraOfVCF == 0 && input$sigAttributionOfVCF == 0){
        return(FALSE)
      }
      return(TRUE)
    }
    
    observeEvent(
      CheckArgumentsForVCF(),
      {
        output$selectSampleForAttribution <- renderUI(
          { 
            sample.names <- colnames(list.of.catalogs[[1]])
            selectInput(inputId = "selectedSampleForAttribution",
                        label = "Select sample",
                        choices = sample.names)
          }
        )
        
        output$selectCancerType <- renderUI(
          {
            cancer.types <-
              c("Unknown", colnames(CancerTypeToExposureStatData()))
            selectInput(inputId = "selectedCancerType",
                        label = "Select cancer type",
                        choices = cancer.types)
          }
        )
        
        output$chooseCatalogType <- renderUI(
          {
            selectizeInput(inputId = "selectCatalogType", 
                        label = "Select catalog type",
                        choices = c("SBS96", "SBS192", "DBS78", "ID"),
                        options = list(
                          placeholder = 'Please select an option below',
                          onInitialize = I('function() { this.setValue(""); }')
                        )
                        )
          }
        )
        
        shinyjs::show(id = "chooseCatalogType")
        
        # Hide the UI element previously generated by uploading spectra
        output$uploadedCatalogType <- NULL
        shinyjs::hide(id = "uploadedCatalogType")
        
        # Hide the previous signature aetiology table
        output$sigAetiologyTable <- NULL
        shinyjs::hide(id = "sigAetiologyTable")
      }) # end of observeEvent
    
    
    observeEvent(input$selectedCancerType, {
      req(input$selectedCancerType)
      
      # Don't update the input.catalog.type using input$selectCatalogType
      # We are doing analysis using uploaded spectra, however input$selectCatalogType
      # is related to uploaded VCF
      if (is.null(catalog)) {
        catalog.name <- paste0("cat", input.catalog.type())
        catalog <<- list.of.catalogs[[catalog.name]]
      }
      ShowPreselectedSigs(input, output, input.catalog.type())
      
      # Show the actionButton for user to add more signatures
      output$addSig <- renderUI(
        actionButton(inputId = "addMoreSigs", label = "Add more signatures",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      )
      
    }, ignoreInit = TRUE) # end of observeEvent
    
    observeEvent(input$selectCatalogType, {
      req(input$selectCatalogType)
      
      input.catalog.type(input$selectCatalogType) 
      if (is.null(catalog)) {
        catalog.name <- paste0("cat", input.catalog.type())
        catalog <<- list.of.catalogs[[catalog.name]]
      }
      ShowPreselectedSigs(input, output, input.catalog.type())
      
      # Show the actionButton for user to add more signatures
      output$addSig <- renderUI(
        actionButton(inputId = "addMoreSigs", label = "Add more signatures",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      )
      
    }, ignoreInit = TRUE) # end of observeEvent
    
    # Update the catalog used for attribution if user selects another catalog type
    # for analysis
    observeEvent(input$selectCatalogType, {
      catalog.name <- paste0("cat", input$selectCatalogType)
      catalog <<- list.of.catalogs[[catalog.name]]
    })
    
    observeEvent(input$addMoreSigs, {
      output$chooseMoreSigs <- renderUI(
        {
          shinyWidgets::pickerInput(inputId = "selectedMoreSigs",
                                    label = "Choose more signatures",  
                                    choices = choose.more.sigs(),
                                    options = shinyWidgets::pickerOptions(
                                      actionsBox = TRUE,
                                      dropupAuto = FALSE,
                                      `live-search`=TRUE
                                    ), 
                                    multiple = TRUE)
        }
      )
    })
    
    # Create a reactive expression to determine the signatures used for
    # attribution and then show the signatures aetiology information table
    sigsForAttribution <- reactive({
      list(input$selectCatalogType ,input$preselectedSigs, input$selectedMoreSigs)
      return(c(input$preselectedSigs, input$selectedMoreSigs))
    })
    
    observeEvent(sigsForAttribution(), {
      
      if (is.null(input.catalog.type())) {
        return()
      }
      
      if (input.catalog.type() == "SBS96") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.SBS96.sig.links), sigsForAttribution())
      } else if (input.catalog.type() == "SBS192") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.SBS192.sig.links), sigsForAttribution())
      } else if (input.catalog.type() == "DBS78") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.DBS78.sig.links), sigsForAttribution())
      } else if (input.catalog.type() == "ID") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.ID.sig.links), sigsForAttribution())
      } else {
        sigs.in.correct.order <- NULL
      }
      
      
      dat <<- PrepareSigsAetiologyTable(input.catalog.type(), 
                                        input.ref.genome(), input.region())
      dat1 <- dat[sigs.in.correct.order, ]
      plotdata$dat <<- dat1
      
      df <- PrepareThumbnailForSample(input = input, catalog = catalog,
                                      input.catalog.type = input.catalog.type())
      
      dat2 <- rbind(df, dat1)
      output$sigAetiologyTable <- DT::renderDataTable({
        DT::datatable(dat2, escape = FALSE, rownames = FALSE,
                      colnames = c("Name", "Spectrum", "Proposed etiology"), 
                      options = list(lengthMenu = c(25, 50, 75), 
                                     pageLength = 25,
                                     language = list(
                                       search = "Search in signatures ID and etiologies:"
                                     )
                      )
        )
      })
      
      # If no signatures were selected, hide the signature etiology table
      if (is.null(sigsForAttribution())) {
        output$sigAetiologyTable <- NULL
        shinyjs::hide(id = "sigAetiologyTable")
      } else {
        shinyjs::show(id = "sigAetiologyTable")
      }
      
      # Don't show "Analyze" button if no signatures were selected in the beginning
      if (!is.null(sigsForAttribution())) {
        output$analysisButton <- renderUI({
          actionButton(inputId = "startAnalysis", label = "Analyze",
                       style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
        })
      } # Need to use ignoreNULL here otherwise sigAetiologyTable will not change
        # if user deselect all signatures
    }, ignoreNULL = FALSE,  ignoreInit = TRUE)
    
    # When user changes input for sample selected and catalog type
    # hide the previous signatuer aetiology table
    hideAetiologyTable <- reactive(
      list(input$selectedSampleForAttribution,
           input$selectCatalogType)
    )
    
    # Asynchronous programming starts from here
    observeEvent(input$startAnalysis, {
      if(length(sigsForAttribution()) == 0) {
        showNotification(ui = "Error:", 
                         action = "No signatures selected for attribution analysis",
                         type = "error", duration = NULL)
        return()
      }
      
      # Don't do anything if in the middle of a run
      if(running())
        return(NULL)
      
      running(TRUE)
      spect <- catalog[, input$selectedSampleForAttribution, drop = FALSE]
      catalog.type <- input.catalog.type()
      cancer.type <- input$selectedCancerType
      
      if (catalog.type == "ID") {
        sig.universe <<- 
          COSMIC.v3.hg19.genome.ID.sigs[, sigsForAttribution(), drop = FALSE]
      } else {
        
        sig.universe <<- total.signatures()[, sigsForAttribution(), drop = FALSE]
      }
      
      # Do the first round of cut-off if there are many signatures in the beginning
      if (ncol(sig.universe) > 15) {
        retval <- 
          mSigAct::OptimizeExposureQPBootstrap(spectrum = spect,
                                               signatures = sig.universe, 
                                               mc.cores = AdjustNumberOfCores(50))
        sig.universe <<- sig.universe[, names(retval$exposure), drop = FALSE]
      }
      
      sigs.prop <- mSigAct::ExposureProportions(
        mutation.type = catalog.type,
        cancer.type   = cancer.type,
        all.sigs      = sig.universe,
        must.include  = colnames(sig.universe))
      
      # Create a Progress object
      progress <- ipc::AsyncProgress$new(session, min = 0, max = 1,
                                         message = "Analysis in progress",
                                         detail = "This may take a while...")
      
      result_val(NULL)
      #future.planned <<- TRUE
      
      # future::multicore is not supported on Windows
      if (.Platform$OS.type != "windows") {
        if (!future.planned) {
          future.planning.time <- system.time(
            future::plan(future::multicore(workers = min(64, future::availableCores())))
          )
          print(future.planning.time)
        }
        
        future.planned <<- TRUE
      }
      
      fut <- future::future(
        {
          # Close the progress when this reactive exits (even if there's an error)
          # on.exit(progress$close())
          
          # Create a callback function to update progress. Each time this is called, it
          # will increase the progress by that value and update the detail
          updateProgress <- function(value = NULL, detail = NULL) {
            progress$inc(amount = value, detail = detail)
            interruptor$execInterrupts()
          }
          
          # In some cases the first phase of the search, which is a 
          # global search, might find different optima depending
          # on the random seed, so the purpose of set.seed with L'Ecuyer-CMRG
          # is to get reproducible results across mclapply / mcparallel 
          # (https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mcparallel.html).
          set.seed(102119, kind = "L'Ecuyer-CMRG")
          
          retval <- mSigAct::MAPAssignActivity1(
            spect = spect,
            sigs = sig.universe,
            sigs.presence.prop = sigs.prop,
            max.level = length(sigs.prop) - 1,
            p.thresh = 0.01,
            m.opts = mSigAct::DefaultManyOpts(),
            max.mc.cores = AdjustNumberOfCores(50),
            progress.monitor = updateProgress
          )
          
          return(retval)
        }, seed = TRUE) %...>% {
          retval <- .
          
          plotdata$retval <<- retval
          
          if (retval$success == FALSE || is.null(retval$success)) {
            showNotification(ui = "Message:", 
                             action = retval$messages,
                             type = "message", duration = NULL)
            return()
          } else {
            MAP.best.exp <- retval$MAP
            
            reconstructed.catalog0 <- retval$MAP.recon
            
            cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
            
            colnames(reconstructed.catalog0) <- 
              paste0("reconstructed (cosine similarity = ", cossim, ")")
            
            plotdata$cossim <<- cossim
            plotdata$spect <<- spect
            plotdata$reconstructed.catalog <<- reconstructed.catalog0
            plotdata$sig.universe <<- sig.universe
            plotdata$best.MAP.exp <<- MAP.best.exp
            
            retval <- PrepareAttributionResults2(input, output, session, 
                                                 input.catalog.type(), 
                                                 plotdata)
            attribution.results <<- retval$attribution.results
          }
        } %...>% result_val
      
      # Show notification on error or user interrupt
      fut <- promises::catch(fut,
                             function(e){
                               result_val(NULL)
                               print(e$message)
                               showNotification(e$message, 
                                                duration = NULL)
                             })
      
      # When done with analysis, remove progress bar
      fut <- promises::finally(fut, function(){
        progress$close()
        running(FALSE) # Declare done with run
      })
      
      # Return something other than the future so we don't block the UI
      NULL
    })
    
    # Send interrupt signal to future
    observeEvent(input$cancel,{
      if(running()) {
        interruptor$interrupt("Task cancelled")
      }
    })
    
  }) # end of tryCatch
}
  
