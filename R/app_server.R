# Cannot use plan(multicore), otherwise the progress bar for asynchronous
# process will not work properly
future::plan(future::multisession(workers = min(64, future::availableCores())))

#' @import mSigAct
#' @import promises
#' @import ipc
#' @import shiny
#' @import shinydashboard
app_server <- function(input, output, session) {
  # List the first level callModules here
  tryCatch({
    fut <- NULL
    result_val <- reactiveVal()
    running <- reactiveVal(FALSE)
    interruptor <- ipc::AsyncInterruptor$new()
    
    # Create reactiveValues object
    # and set flag to 0 to prevent errors with adding NULL
    rv <- reactiveValues(downloadFlag = 0)
    
    # Increase the file upload limit to 100MB
    options(shiny.maxRequestSize=100*1024^2)
    
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
    
    # Create a variable which can be used to store the uploaded catalog later
    catalog <- NA
    
    list.of.catalogs <- NA
    
    choose.more.sigs <- catalog.path <- input.catalog.type <- NA
    
    showSBS192Catalog <- TRUE
    
    attribution.results <- FALSE
    
    plotdata <- list(cossim = NULL, spect = NULL, 
                     reconstructed.catalog = NULL,
                     sig.universe = NULL, 
                     QP.best.MAP.exp = NULL,
                     dat = NULL)
    
    dat <- data.frame(name = character(), spectrum = character(), 
                      proposed.aetiology = character())
    
    #######################################################################
    # Functions related to UploadVCFUI(), the first tab
    #######################################################################
    # When user uploads VCF files, then show the action button to ge
    
    
    
    # Download sample VCFs when user clicks the button
    output$downloadsampleVCFs <- downloadHandler(
      filename = function() {
        "mSigAct-example-VCFs.zip"
      },
      content = function(file) {
        PrepareExampleVCFs(file)
      })
    
    # Download sample catalogs when user clicks the button
    output$downloadSampleSpectra <- downloadHandler(
      filename = function() {
        "mSigAct-example-spectra.zip"
      },
      content = function(file) {
        PrepareExampleSpectra(file)
      })
    
    # Run analysis on sample Strelka SBS VCFs when user clicks the button
    output$runstrelkasbsvcfs <- downloadHandler(
      filename = function() {
        "mSigAct-test-run-Strelka-VCFs-output.zip"
      },
      content = function(file) {
        results <- RunICAMSOnSampleStrelkaVCFs(session, output, file, ids)
        list.of.catalogs <<- results$counts
        # When the downloadHandler function runs, increment rv$downloadFlag
        rv$downloadFlag <- rv$downloadFlag + 1
      })
    
    # Run analysis on sample Mutect VCFs when user clicks the button
    output$runmutectvcfs <- downloadHandler(
      filename = function() {
        "mSigAct-test-run-Mutect-VCFs-output.zip"
      },
      content = function(file) {
        results <- RunICAMSOnSampleMutectVCFs(session, output, file, ids)
        list.of.catalogs <<- results$counts
        # When the downloadHandler function runs, increment rv$downloadFlag
        rv$downloadFlag <- rv$downloadFlag + 1
      })
    
    if (FALSE) {
      output$SBS96plot <- NULL
      output$SBS192plot <- NULL
      output$SBS1536plot <- NULL
      output$DBS78plot <- NULL
      output$DBS136plot <- NULL
      output$DBS144plot <- NULL
      output$IDplot <- NULL
    }
    
    # When user submits VCF to analysis, then the program will try to
    # generate a zip archive based on the input files and parameters
    output$download <- downloadHandler(
      
      filename = function() {
        errors <- CheckInputsForVCF(input)
        ids$error <<- append(ids$error, AddErrorMessage(errors))
        paste0(input$zipfile.name, ".zip")
      },
      
      content = function(file) {
        req(input$variantCaller, input$ref.genome, input$region, input$vcf.files)
        if (input$variantCaller == "unknown") {
          if (is.null(input$mergeSBS)) {
            return()
          }
        }
        
        result <- ProcessVCFs(input, output, file, ids)
        retval <<- result$retval
        old.error.ids <- ids$error
        ids <<- result$ids
        ids$error <- append(ids$error, old.error.ids)
        
        list.of.catalogs <<- retval$counts
        #density.catalog <- retval$density
        
        # When the downloadHandler function runs, increment rv$downloadFlag
        rv$downloadFlag <- rv$downloadFlag + 1

      })
    
    
    observeEvent(rv$downloadFlag, {
      output$showSpectraFromVCF <- renderUI({
        actionButton(inputId = "showSpectraOfVCF", label = "Show spectra",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
      })
      
      output$sigAttributionFromVCF <- renderUI({
        actionButton(inputId = "sigAttributionOfVCF", label = "Signature attribution",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
      })
      
      output$selectSampleFromUploadedVCF <- renderUI(
        {
          sample.names <- colnames(list.of.catalogs[[1]])
          radioButtons(inputId = "sampleNameFromUploadedVCF",
                       label = "Select sample from uploaded VCF",
                       choices = sample.names,
                       selected = character(0))
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
      if (input.catalog.type %in% c("SBS96", "SBS192", "DBS78", "ID")) {
        shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
      }
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
    ########################################################################
    ShowTwoButtons <- function() {
      output$showSpectraFromCatalog <- renderUI(
        actionButton(inputId = "showSpectraFromCatalog", label = "Show spectra",
                     style="color: #fff; background-color: #337ab7;
                          border-color: #2e6da4;")
      )
      output$sigAttributionFromCatalog <- renderUI(
        actionButton(inputId = "sigAttributionFromCatalog", 
                     label = "Signature attribution",
                     style="color: #fff; background-color: #337ab7;
                          border-color: #2e6da4;"))
    }
    
    observeEvent(input$preloadSBS96Spectra, {
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "genome")
      catalog.path <<- system.file("extdata/SBS96-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type <<- "SBS96"
      
      ShowTwoButtons()
      
    })
    
    observeEvent(input$preloadSBS192Spectra, {
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "transcript")
      catalog.path <<- system.file("extdata/SBS192-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type <<- "SBS192"
      
      ShowTwoButtons()
      
    })
    
    observeEvent(input$preloadDBS78Spectra, {
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "genome")
      catalog.path <<- system.file("extdata/DBS78-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type <<- "DBS78"
      
      ShowTwoButtons()
      
    })
    
    observeEvent(input$preloadIDSpectra, {
      
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "ref.genome2",
                                      selected = "hg19")
      shinyWidgets::updatePickerInput(session = session,
                                      inputId = "region2",
                                      selected = "genome")
      catalog.path <<- system.file("extdata/ID-mSigAct-example-spectra.csv", 
                                   package = "mSigAct.server")
      input.catalog.type <<- "ID"
      
      ShowTwoButtons()
      
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
    })
    
    
    observeEvent(CheckArgumentsForSpectra(), {
      
      if (TwoActionButtonsClicked(input) == FALSE) {
        return()
      } else {
        output$removeButton2 <- renderUI(
          actionButton(inputId = "remove2",
                       label = "Remove notifications"))
      }
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
    observeEvent(input$showSpectraFromCatalog, {
      
      if (TwoActionButtonsClicked(input) == FALSE) {
        return()
      } else if (showSBS192Catalog == FALSE) {
        return()
      } else {
        req(input$ref.genome2, input$region2)
        # Hide the widgets for previous uploaded VCF
        output$spectraPlotFromVCF <- NULL
        output$selectSampleFromUploadedVCF <- NULL
        shinyjs::hide(id = "spectraPlotFromVCF")
        shinyjs::hide(id = "selectSampleFromUploadedVCF")
        
        output$selectSampleFromVCFForAttribution <- NULL
        shinyjs::hide(id = "selectSampleFromVCFForAttribution")
        
        output$selectCancerTypeOfVCF <- NULL
        shinyjs::hide(id = "selectCancerTypeOfVCF")
        
        output$chooseCatalogType <- NULL
        shinyjs::hide(id = "chooseCatalogType")
        
        output$chooseSigSubsetForVCF <- NULL
        shinyjs::hide(id = "chooseSigSubsetForVCF")
        
        output$addSigForVCF <- NULL
        shinyjs::hide(id = "addSigForVCF")
        
        output$chooseMoreSigsForVCF <- NULL
        shinyjs::hide(id = "chooseMoreSigsForVCF")
        
        output$analysisButtonForVCF <- NULL
        shinyjs::hide(id = "analysisButtonForVCF")
        
        output$mytableForVCF <- NULL
        shinyjs::hide(id = "mytableForVCF")
        
        shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
        
        shinyjs::show(id = "mytable")
        if (input.catalog.type %in% c("SBS96", "SBS192", "DBS78", "ID")) {
          shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
        }
        shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                       selected = "showSpectraTab")
      }
    })
    
    observeEvent(input$sigAttributionFromCatalog, {
      
      if (TwoActionButtonsClicked(input) == FALSE) {
        return()
      } else if (showSBS192Catalog == FALSE) {
        return()
      } else {
        req(input$ref.genome2, input$region2)
        if (!input.catalog.type %in% c("SBS96", "SBS192", "DBS78", "ID")) {
          showNotification(ui = "Error:", 
                           action = paste0("Can only do signature attribution ", 
                                           "for SBS96, SBS192, DBS78 and ID"),
                           type = "error")
          return()
        }
        # Hide the widgets for previous uploaded VCF
        shinyjs::hide(id = "spectraPlotFromVCF")
        shinyjs::hide(id = "selectSampleFromUploadedVCF")
        
        shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
        shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
        shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                       selected = "sigAttributionTab2")
      }
    })
    
    observeEvent(CheckArgumentsForSpectra(), {
      
      req(input$ref.genome2, input$region2)
      
      # Delete the previous spectra plot
      output$spectraPlotFromCatalog <- NULL
      
      
      if (showSBS192Catalog == FALSE) {
        return()
      }
      
      output$selectSampleFromUploadedCatalog <-
        renderUI(
          {
            catalog <<- ICAMS::ReadCatalog(file = catalog.path,
                                           ref.genome = input$ref.genome2,
                                           region = input$region2)
            
            if (nrow(catalog) == 96) {
              input.catalog.type <<- "SBS96"
            } else if (nrow(catalog) == 192) {
              input.catalog.type <<- "SBS192"
            } else if (nrow(catalog) == 1536) {
              input.catalog.type <<- "SBS1536"
            } else if (nrow(catalog) == 78) {
              input.catalog.type <<- "DBS78"
            } else if (nrow(catalog) == 136) {
              input.catalog.type <<- "DBS136"
            } else if (nrow(catalog) == 144) {
              input.catalog.type <<- "DBS144"
            } else if (nrow(catalog) == 83) {
              input.catalog.type <<- "ID"
            }
            
            sample.names <- colnames(catalog)
            
            tagList(
              radioButtons(inputId = "selectedSampleFromUploadedCatalog",
                           label = "Select spectrum to view",
                           choices = sample.names,
                           selected = character(0)),
              
              if (input.catalog.type %in% c("SBS96", "SBS192", "DBS78", "ID")) {
                actionButton(inputId = "clickToSigAttribution",
                             label = "Signature attribution",
                             style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
              }
              
            )
          }
        )
    })
    
    # When user clicks the action button on Show spectra page, direct user to the relevant tab
    observeEvent(input$clickToSigAttribution, {
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab2")
    })
    
    # When user submit new catalog for analysis, remove the previous plots
    observeEvent(input$showSpectraFromCatalog, {
      output$SBS96plot <- NULL
      output$SBS192plot <- NULL
      output$SBS1536plot <- NULL
      output$DBS78plot <- NULL
      output$DBS136plot <- NULL
      output$DBS144plot <- NULL
      output$IDplot <- NULL
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
      
      
      if (input.catalog.type == "SBS96") {
        height <- 230
        width <- 800
      } else if (input.catalog.type == "SBS192") {
        height <- 250
        width <- 800
      } else if (input.catalog.type == "SBS1536") {
        height <- 800
        width <- 800
      } else if (input.catalog.type == "DBS78") {
        height <- 250
        width <- 800
      } else if (input.catalog.type == "DBS136") {
        height <- 500
        width <- 700
      } else if (input.catalog.type == "DBS144") {
        height <- 350
        width <- 350
      } else if (input.catalog.type == "ID") {
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
    
    CheckArgumentsForSpectra <- reactive({
      list(input$showSpectraFromCatalog, input$sigAttributionFromCatalog)
    })
    
    TwoActionButtonsClicked <- function(input) {
      if(is.null(input$showSpectraFromCatalog) && is.null(input$sigAttributionFromCatalog)){
        return(FALSE)
      }
      if(input$showSpectraFromCatalog == 0 && input$sigAttributionFromCatalog == 0){
        return(FALSE)
      }
      return(TRUE)
    }
    
    observeEvent(
      CheckArgumentsForSpectra(),
      {
        req(input$ref.genome2, input$region2)
        
        if (showSBS192Catalog == FALSE) {
          return()
        }
        
        
        output$selectSampleFromCatalogForAttribution2 <- renderUI(
          { 
            catalog <<- ICAMS::ReadCatalog(file = catalog.path,
                                           ref.genome = input$ref.genome2,
                                           region = input$region2)
            
            if (nrow(catalog) == 96) {
              input.catalog.type <<- "SBS96"
            } else if (nrow(catalog) == 192) {
              input.catalog.type <<- "SBS192"
            } else if (nrow(catalog) == 1536) {
              input.catalog.type <<- "SBS1536"
            } else if (nrow(catalog) == 78) {
              input.catalog.type <<- "DBS78"
            } else if (nrow(catalog) == 136) {
              input.catalog.type <<- "DBS136"
            } else if (nrow(catalog) == 144) {
              input.catalog.type <<- "DBS144"
            } else if (nrow(catalog) == 83) {
              input.catalog.type <<- "ID"
            }
            
            sample.names <- colnames(catalog)
            selectInput(inputId = "selectedSampleFromCatalogForAttribution2",
                        label = "Select sample from uploaded spectra",
                        choices = sample.names)
          }
        )
        
        output$selectCancerType2 <- renderUI(
          {
            cancer.types <-
              c("Unknown", colnames(CancerTypeToExposureStatData()))
            selectInput(inputId = "selectedCancerType2",
                        label = "Select cancer type",
                        choices = cancer.types)
          }
        )
        
        
        output$uploadedCatalogType <- renderUI(
          {
            p(tags$b("Mutation type: "), input.catalog.type)
          }
        )
      })
    
    observeEvent(input$selectedCancerType2, {
      if (!input.catalog.type %in% c("SBS96", "SBS192", "DBS78", "ID")) {
        showNotification(ui = "Error:", 
                         action = paste0("Can only do signature attribution ", 
                                         "for SBS96, SBS192, DBS78 and ID"),
                         type = "error")
        return()
      }
      
      if (input$selectedCancerType2 != "Unknown") {
        output$chooseSigSubsetForSampleFromCatalog2 <- renderUI(
          { 
            
            if (input.catalog.type == "SBS96") {
              sig.universe <- colnames(COSMIC.v3.genome.SBS96.sigs)
            } else {
              sig.universe <- 
                colnames(PCAWG7::signature[["genome"]][[input.catalog.type]])
            }
            
            if (input$selectedCancerType2 == "Unknown") {
              selected.sig.universe <- NULL
            } else {
              tmp <- CancerTypeToSigSubset(cancer.type = input$selectedCancerType2,
                                           tumor.cohort = "PCAWG",
                                           sig.type = input.catalog.type,
                                           region = "genome")
              selected.sig.universe0 <- colnames(tmp)
              
              # Exclude possible artifact signatures
              possible.artifacts <- mSigAct::PossibleArtifacts()
              
              selected.sig.universe1 <- 
                setdiff(selected.sig.universe0, possible.artifacts)
              
              # Exclude rare signatures
              rare.sigs <- mSigAct::RareSignatures()
              selected.sig.universe <-
                setdiff(selected.sig.universe1, rare.sigs)
            }
            
            choose.more.sigs <<- setdiff(sig.universe, selected.sig.universe)
            shinyWidgets::pickerInput (inputId = "preselectedSigs",
                                       label = paste0("These signatures were preselected based ",  
                                                      "on cancer type"),
                                       choices = selected.sig.universe,
                                       selected = selected.sig.universe,
                                       options = shinyWidgets::pickerOptions(
                                         actionsBox = TRUE,
                                         dropupAuto = FALSE
                                       ), 
                                       multiple = TRUE
            )
            
            
            
            
          }
        )
        output$addSig <- renderUI(
          actionButton(inputId = "addMoreSigs", label = "Add more signatures",
                       style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
        )
        
        if (input.catalog.type == "SBS96") {
          dat <<- data.frame(
            name = paste0("<a href='", COSMIC.v3.SBS.sig.links, "' target='_blank'>", 
                          rownames(COSMIC.v3.SBS.sig.links),  "</a>"), 
            spectrum = paste0('<img src="SBS96/', rownames(COSMIC.v3.SBS.sig.links), '.PNG"',
                              ' height="52"></img>'),
            proposed.aetiology = SBS.aetiology)
        } else if (input.catalog.type == "SBS192") {
          dat <<- data.frame(
            name = paste0("<a href='", COSMIC.v3.SBS.sig.links, "' target='_blank'>", 
                          rownames(COSMIC.v3.SBS.sig.links),  "</a>"), 
            spectrum = paste0('<img src="SBS192/', rownames(COSMIC.v3.SBS.sig.links), '.PNG"',
                              ' height="52"></img>'),
            proposed.aetiology = SBS.aetiology)
        } else if (input.catalog.type == "DBS78") {
          dat <<- data.frame(
            name = paste0("<a href='", COSMIC.v3.DBS.sig.links, "' target='_blank'>", 
                          rownames(COSMIC.v3.DBS.sig.links),  "</a>"), 
            spectrum = paste0('<img src="DBS78/', rownames(COSMIC.v3.DBS.sig.links), '.PNG"',
                              ' height="52"></img>'),
            proposed.aetiology = DBS.aetiology)
        } else if (input.catalog.type == "ID") {
          dat <<- data.frame(
            name = paste0("<a href='", COSMIC.v3.ID.sig.links, "' target='_blank'>", 
                          rownames(COSMIC.v3.ID.sig.links),  "</a>"), 
            spectrum = paste0('<img src="ID/', rownames(COSMIC.v3.ID.sig.links), '.PNG"',
                              ' height="52"></img>'),
            proposed.aetiology = ID.aetiology)
        } 
        
        #rownames(dat) <- names(COSMIC.v3.SBS.sig.links)
        
        output$mytable <- DT::renderDataTable({
          DT::datatable(dat[input$preselectedSigs, ], escape = FALSE, rownames = FALSE) # HERE)
        })
      }
      
    })
    
    observeEvent(input$addMoreSigs, {
      output$chooseMoreSigs <- renderUI(
        {
          shinyWidgets::pickerInput(inputId = "selectedMoreSigs",
                                    label = "Choose more signatures",  
                                    choices = choose.more.sigs,
                                    options = shinyWidgets::pickerOptions(
                                      actionsBox = TRUE,
                                      dropupAuto = FALSE
                                    ), 
                                    multiple = TRUE)
        }
      )
    })
    
    sigsForAttribution <- reactive({
      sigs.to.show <- c(input$preselectedSigs, input$selectedMoreSigs)
      if (is.na(input.catalog.type)) {
        return()
      }
      
      if (input.catalog.type %in% c("SBS96", "SBS192")) {
        sigs.in.correct.order <- intersect(rownames(COSMIC.v3.SBS.sig.links), sigs.to.show)
      } else if (input.catalog.type == "DBS78") {
        sigs.in.correct.order <- intersect(rownames(COSMIC.v3.DBS.sig.links), sigs.to.show)
      } else if (input.catalog.type == "ID") {
        sigs.in.correct.order <- intersect(rownames(COSMIC.v3.ID.sig.links), sigs.to.show)
      } else {
        return()
      }
      
      return(sigs.in.correct.order)
    })
    
    observeEvent(sigsForAttribution(), {
      
      plotdata$dat <<- dat[sigsForAttribution(), ]
      
      output$mytable <- DT::renderDataTable({
        DT::datatable(dat[sigsForAttribution(), ], escape = FALSE, rownames = FALSE,
                      options = list(lengthMenu = c(25, 50, 75), 
                                     pageLength = 25)) 
      })
      
      output$analysisButton <- renderUI({
        actionButton(inputId = "startAnalysis", label = "Analyze",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
      })
      
    })
    
    observeEvent(input$startAnalysis, {
      
      if(length(sigsForAttribution()) == 0) {
        showNotification(ui = "Error:", 
                         action = "No signatures selected for attribution analysis",
                         type = "error")
        return()
      }
      
      #Don't do anything if in the middle of a run
      if(running())
        return(NULL)
      running(TRUE)
      
      spect <- catalog[, input$selectedSampleFromCatalogForAttribution2, drop = FALSE]
      catalog.type <- input.catalog.type
      cancer.type <- input$selectedCancerType2
      region <- input$region2
      
      if (catalog.type == "SBS192") {
        sig.universe <- 
          PCAWG7::signature[["genome"]][[catalog.type]][, sigsForAttribution(), drop = FALSE]
      } else if (catalog.type == "SBS96") {
        sig.universe <- 
          COSMIC.v3.genome.SBS96.sigs[, sigsForAttribution(), drop = FALSE]
      } else {
        sig.universe <- 
          PCAWG7::signature[[region]][[catalog.type]][, sigsForAttribution(), drop = FALSE]
      }
      
      sigs.prop <- mSigAct::ExposureProportions(mutation.type = catalog.type,
                                                cancer.type = cancer.type,
                                                all.sigs = sig.universe)
      
      if (FALSE) {
        QP.exposure <- 
          GetExposureWithConfidence(catalog = spect,
                                    sig.universe = sig.universe,
                                    num.of.bootstrap.replicates = 10000,
                                    method = decomposeQP,
                                    conf.int = 0.95)
        updated.sig.universe <- sig.universe[ , rownames(QP.exposure)]
        updated.sigs.prop <- sigs.prop[rownames(QP.exposure)]
      }
      
      # Create a Progress object
      progress <- ipc::AsyncProgress$new(session, min = 0, max = 1,
                                         message = "Analysis in progress",
                                         detail = "This may take a while...")
      result_val(NULL)
      
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
          
          AdjustNumberOfCores <- 
            getFromNamespace(x = "Adj.mc.cores", ns = "mSigAct")
          
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
          
          MAP.best.exp <- retval$MAP
          
          QP.exp <- 
            mSigAct::OptimizeExposureQP(spect, 
                                        sig.universe[ , MAP.best.exp$sig.id, 
                                                      drop = FALSE])
          QP.best.MAP.exp <-
            dplyr::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)
          
          r.qp <- mSigAct::ReconstructSpectrum(sig.universe, exp = QP.exp, 
                                               use.sig.names = TRUE)
          reconstructed.catalog0 <- 
            as.catalog(r.qp, ref.genome = input$ref.genome2, 
                       region = input$region2)
          
          cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
          
          colnames(reconstructed.catalog0) <- 
            paste0("reconstructed (cosine similarity = ", cossim, ")")
          reconstructed.catalog <- round(reconstructed.catalog0)
          
          plotdata$cossim <<- cossim
          plotdata$spect <<- spect
          plotdata$reconstructed.catalog <<- reconstructed.catalog
          plotdata$sig.universe <<- sig.universe
          plotdata$QP.best.MAP.exp <<- QP.best.MAP.exp
          
          retval <- PrepareAttributionResults2(input, output, session, 
                                               input.catalog.type, 
                                               plotdata)
          attribution.results <<- retval$attribution.results
        } %...>% result_val
      
      # Show notification on error or user interrupt
      fut <- promises::catch(fut,
                             function(e){
                               result_val(NULL)
                               print(e$message)
                               showNotification(e$message)
                             })
      
      # When done with analysis, remove progress bar
      fut <- promises::finally(fut, function(){
        progress$close()
        running(FALSE) # Declare done with run
      })
      
      # Return something other than the future so we don't block the UI
      NULL
    })
    #######################################################################
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
        #req(input$ref.genome2, input$region2)
        
        #if (showSBS192Catalog == FALSE) {
        #  return()
        #}
        
        output$selectSampleFromVCFForAttribution <- renderUI(
          { 
            sample.names <- colnames(list.of.catalogs[[1]])
            selectInput(inputId = "selectedSampleFromVCFForAttribution",
                        label = "Select sample from uploaded VCF",
                        choices = sample.names)
          }
        )
        
        output$selectCancerTypeOfVCF <- renderUI(
          {
            cancer.types <-
              c("Unknown", colnames(CancerTypeToExposureStatData()))
            selectInput(inputId = "selectedCancerTypeOfVCF",
                        label = "Select cancer type",
                        choices = cancer.types)
          }
        )
        
        
        output$chooseCatalogType <- renderUI(
          {
            selectInput(inputId = "selectCatalogTypeOfVCF", 
                         label = "Select catalog type",
                         choices = c("SBS96", "SBS192", "DBS78", "ID"))
          }
        )
      })
    
    CheckArgumentsForVCFAttribution <- reactive({
      list(input$selectedCancerTypeOfVCF, input$selectCatalogTypeOfVCF)
    })
    
    observeEvent(CheckArgumentsForVCFAttribution(), {
      req(input$selectedCancerTypeOfVCF, input$selectCatalogTypeOfVCF)
      
      input.catalog.type <<- input$selectCatalogTypeOfVCF
      output$chooseSigSubsetForVCF <- renderUI(
        { 
          
          if (input.catalog.type == "SBS96") {
            sig.universe <- colnames(COSMIC.v3.genome.SBS96.sigs)
          } else {
            sig.universe <- 
              colnames(PCAWG7::signature[["genome"]][[input.catalog.type]])
          }
          
          if (input$selectedCancerTypeOfVCF == "Unknown") {
            selected.sig.universe <- NULL
          } else {
            tmp <- CancerTypeToSigSubset(cancer.type = input$selectedCancerTypeOfVCF,
                                         tumor.cohort = "PCAWG",
                                         sig.type = input.catalog.type,
                                         region = "genome")
            selected.sig.universe0 <- colnames(tmp)
            
            # Exclude possible artifact signatures
            possible.artifacts <- mSigAct::PossibleArtifacts()
            
            selected.sig.universe1 <- 
              setdiff(selected.sig.universe0, possible.artifacts)
            
            # Exclude rare signatures
            rare.sigs <- mSigAct::RareSignatures()
            selected.sig.universe <-
              setdiff(selected.sig.universe1, rare.sigs)
          }
          
          choose.more.sigs <<- setdiff(sig.universe, selected.sig.universe)
          shinyWidgets::pickerInput (inputId = "preselectedSigsForVCF",
                                     label = paste0("These signatures were preselected based ",  
                                                    "on cancer type"),
                                     choices = selected.sig.universe,
                                     selected = selected.sig.universe,
                                     options = shinyWidgets::pickerOptions(
                                       actionsBox = TRUE,
                                       dropupAuto = FALSE
                                     ), 
                                     multiple = TRUE
          )
        }
      )
    }
    )
    
    observeEvent(input$preselectedSigsForVCF, {
      output$addSigForVCF <- renderUI(
        actionButton(inputId = "addMoreSigsForVCF", label = "Add more signatures",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
      )
      #if (is.na(req(input$selectCatalogTypeOfVCF))) {
      #  return()
      #}
      if (input.catalog.type == "SBS96") {
        dat <<- data.frame(
          name = paste0("<a href='", COSMIC.v3.SBS.sig.links, "' target='_blank'>", 
                        rownames(COSMIC.v3.SBS.sig.links),  "</a>"), 
          spectrum = paste0('<img src="SBS96/', rownames(COSMIC.v3.SBS.sig.links), '.PNG"',
                            ' height="52"></img>'),
          proposed.aetiology = SBS.aetiology)
      } else if (input.catalog.type == "SBS192") {
        dat <<- data.frame(
          name = paste0("<a href='", COSMIC.v3.SBS.sig.links, "' target='_blank'>", 
                        rownames(COSMIC.v3.SBS.sig.links),  "</a>"), 
          spectrum = paste0('<img src="SBS192/', rownames(COSMIC.v3.SBS.sig.links), '.PNG"',
                            ' height="52"></img>'),
          proposed.aetiology = SBS.aetiology)
      } else if (input.catalog.type == "DBS78") {
        dat <<- data.frame(
          name = paste0("<a href='", COSMIC.v3.DBS.sig.links, "' target='_blank'>", 
                        rownames(COSMIC.v3.DBS.sig.links),  "</a>"), 
          spectrum = paste0('<img src="DBS78/', rownames(COSMIC.v3.DBS.sig.links), '.PNG"',
                            ' height="52"></img>'),
          proposed.aetiology = DBS.aetiology)
      } else if (input.catalog.type == "ID") {
        dat <<- data.frame(
          name = paste0("<a href='", COSMIC.v3.ID.sig.links, "' target='_blank'>", 
                        rownames(COSMIC.v3.ID.sig.links),  "</a>"), 
          spectrum = paste0('<img src="ID/', rownames(COSMIC.v3.ID.sig.links), '.PNG"',
                            ' height="52"></img>'),
          proposed.aetiology = ID.aetiology)
      } 
    })
    
    
    output$mytable <- NULL
    shinyjs::hide(id = "mytable")
    #rownames(dat) <- names(COSMIC.v3.SBS.sig.links)
    
    #output$mytableForVCF <- DT::renderDataTable({
    #  DT::datatable(dat[input$preselectedSigsForVCF, ], escape = FALSE, rownames = FALSE) # HERE)
    #})
    
  
  observeEvent(input$addMoreSigsForVCF, {
    output$chooseMoreSigsForVCF <- renderUI(
      {
        shinyWidgets::pickerInput(inputId = "selectedMoreSigsForVCF",
                                  label = "Choose more signatures",  
                                  choices = choose.more.sigs,
                                  options = shinyWidgets::pickerOptions(
                                    actionsBox = TRUE,
                                    dropupAuto = FALSE
                                  ), 
                                  multiple = TRUE)
      }
    )
  })
  
  sigsForAttributionVCF <- reactive({
    sigs.to.show <- c(input$preselectedSigsForVCF, input$selectedMoreSigsForVCF)
    if (is.na(input.catalog.type)) {
      return()
    }
    
    if (input.catalog.type %in% c("SBS96", "SBS192")) {
      sigs.in.correct.order <- intersect(rownames(COSMIC.v3.SBS.sig.links), sigs.to.show)
    } else if (input.catalog.type == "DBS78") {
      sigs.in.correct.order <- intersect(rownames(COSMIC.v3.DBS.sig.links), sigs.to.show)
    } else if (input.catalog.type == "ID") {
      sigs.in.correct.order <- intersect(rownames(COSMIC.v3.ID.sig.links), sigs.to.show)
    } else {
      return()
    }
    
    return(sigs.in.correct.order)
  })
  
  observeEvent(sigsForAttributionVCF(), {
    
    plotdata$dat <<- dat[sigsForAttributionVCF(), ]
    
    output$mytableForVCF <- DT::renderDataTable({
      DT::datatable(dat[sigsForAttributionVCF(), ], escape = FALSE, rownames = FALSE,
                    options = list(lengthMenu = c(25, 50, 75), 
                                   pageLength = 25)) 
    })
    
    output$analysisButtonForVCF <- renderUI({
      actionButton(inputId = "startAnalysisForVCF", label = "Analyze",
                   style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
    })
    
  })
  
  observeEvent(input$startAnalysisForVCF, {
    
    if(length(sigsForAttributionVCF()) == 0) {
      showNotification(ui = "Error:", 
                       action = "No signatures selected for attribution analysis",
                       type = "error")
      return()
    }
    
    #Don't do anything if in the middle of a run
    if(running())
      return(NULL)
    running(TRUE)
    
    catalog.name <- paste0("cat", input$selectCatalogTypeOfVCF)
    catalog <- list.of.catalogs[[catalog.name]]
    
    spect <- catalog[, input$selectedSampleFromVCFForAttribution, drop = FALSE]
    catalog.type <- input.catalog.type
    cancer.type <- input$selectedCancerTypeOfVCF
    region <- input$region
    
    if (catalog.type == "SBS192") {
      sig.universe <- 
        PCAWG7::signature[["genome"]][[catalog.type]][, sigsForAttributionVCF(), drop = FALSE]
    } else if (catalog.type == "SBS96") {
      sig.universe <- 
        COSMIC.v3.genome.SBS96.sigs[, sigsForAttributionVCF(), drop = FALSE]
    } else {
      sig.universe <- 
        PCAWG7::signature[[region]][[catalog.type]][, sigsForAttributionVCF(), drop = FALSE]
    }
    
    sigs.prop <- mSigAct::ExposureProportions(mutation.type = catalog.type,
                                              cancer.type = cancer.type,
                                              all.sigs = sig.universe)
    
    if (FALSE) {
      QP.exposure <- 
        GetExposureWithConfidence(catalog = spect,
                                  sig.universe = sig.universe,
                                  num.of.bootstrap.replicates = 10000,
                                  method = decomposeQP,
                                  conf.int = 0.95)
      updated.sig.universe <- sig.universe[ , rownames(QP.exposure)]
      updated.sigs.prop <- sigs.prop[rownames(QP.exposure)]
    }
    
    # Create a Progress object
    progress <- ipc::AsyncProgress$new(session, min = 0, max = 1,
                                       message = "Analysis in progress",
                                       detail = "This may take a while...")
    result_val(NULL)
    
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
        
        AdjustNumberOfCores <- 
          getFromNamespace(x = "Adj.mc.cores", ns = "mSigAct")
        
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
        
        MAP.best.exp <- retval$MAP
        
        QP.exp <- 
          mSigAct::OptimizeExposureQP(spect, 
                                      sig.universe[ , MAP.best.exp$sig.id, 
                                                    drop = FALSE])
        QP.best.MAP.exp <-
          dplyr::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)
        
        r.qp <- mSigAct::ReconstructSpectrum(sig.universe, exp = QP.exp, 
                                             use.sig.names = TRUE)
        reconstructed.catalog0 <- 
          as.catalog(r.qp, ref.genome = input$ref.genome, 
                     region = input$region)
        
        cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
        
        colnames(reconstructed.catalog0) <- 
          paste0("reconstructed (cosine similarity = ", cossim, ")")
        reconstructed.catalog <- round(reconstructed.catalog0)
        
        plotdata$cossim <<- cossim
        plotdata$spect <<- spect
        plotdata$reconstructed.catalog <<- reconstructed.catalog
        plotdata$sig.universe <<- sig.universe
        plotdata$QP.best.MAP.exp <<- QP.best.MAP.exp
        
        retval <- PrepareAttributionResultsVCF(input, output, session, 
                                             input.catalog.type, 
                                             plotdata)
        attribution.results <<- retval$attribution.results
      } %...>% result_val
    
    # Show notification on error or user interrupt
    fut <- promises::catch(fut,
                           function(e){
                             result_val(NULL)
                             print(e$message)
                             showNotification(e$message)
                           })
    
    # When done with analysis, remove progress bar
    fut <- promises::finally(fut, function(){
      progress$close()
      running(FALSE) # Declare done with run
    })
    
    # Return something other than the future so we don't block the UI
    NULL
  })
    
    
    
    ###########################################################################
    

    
    # When user clicks the "Remove notifications" button, all the previous
    # notifications(error, warning or message) will be removed
    observeEvent(input$remove, {
      RemoveAllNotifications(ids)
    })
    
    observeEvent(input$remove2, {
      RemoveAllNotifications(ids)
    })
  })
}
  
