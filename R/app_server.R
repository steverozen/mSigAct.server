# Cannot use plan(multicore), otherwise the progress bar for asynchronous
# process will not work properly
#future::plan(future::multisession)

#' @import mSigAct
#' @import promises
#' @import ipc
#' @import shiny
#' @import shinydashboard
app_server <- function(input, output, session) {
  # List the first level callModules here
  
  addResourcePath(prefix = "results", directoryPath = tempdir())
  
  fut <- NULL
  result_val <- reactiveVal()
  running <- reactiveVal(FALSE)
  interruptor <- ipc::AsyncInterruptor$new()

  # Increase the file upload limit to 100MB
  options(shiny.maxRequestSize=100*1024^2)
  
  plotdata <- reactiveValues(cossim = NULL, spect = NULL, 
                             reconstructed.catalog = NULL,
                             sig.universe = NULL, QP.best.MAP.exp = NULL)
  
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
  
  choose.more.sigs <- catalog.path <- input.catalog.type <- NA
  
  showSBS192Catalog <- TRUE
  
  attribution.results <- FALSE
  
  dat <- data.frame(name = character(), spectrum = character(), 
                    proposed.aetiology = character())

  # Create reactiveValues object
  # and set flag to 0 to prevent errors with adding NULL
  #rv <- reactiveValues(downloadFlag = 0)

  # Download sample VCFs when user clicks the button
  output$downloadsampleVCFs <- downloadHandler(
    filename = function() {
      "mSigAct-sample-VCFs.zip"
    },
    content = function(file) {
      PrepareSampleVCFs(file)
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
      "mSigAct-test-run-Strelka-SBS-VCFs-output.zip"
    },
    content = function(file) {
      ids <<-
        RunICAMSOnSampleStrelkaSBSVCFs(output, file, ids)
    })

  # Run analysis on sample Mutect VCFs when user clicks the button
  output$runmutectvcfs <- downloadHandler(
    filename = function() {
      "mSigAct-test-run-Mutect-VCFs-output.zip"
    },
    content = function(file) {
      ids <<-
        RunICAMSOnSampleMutectVCFs(output, file, ids)
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
      if (input$vcftype == "strelka.sbs") {
        # Generate a zip archive from Strelka SBS VCFs and
        # update the notification ids for errors, warnings
        # and messages

        result <- ProcessStrelkaSBSVCFs(input, output, file, ids)
        retval <<- result$retval
        old.error.ids <- ids$error
        ids <<- result$ids
        ids$error <- append(ids$error, old.error.ids)

        counts.catalog <- retval$counts
        density.catalog <- retval$density

        if (FALSE) {
          # When the downloadHandler function runs, increment rv$downloadFlag
          rv$downloadFlag <- rv$downloadFlag + 1

          if(rv$downloadFlag > 0){  # trigger event whenever the value of rv$downloadFlag changes
            output$SBS96plot <- NULL
            output$SBS192plot <- NULL
            output$SBS1536plot <- NULL
            output$DBS78plot <- NULL
            output$DBS136plot <- NULL
            output$DBS144plot <- NULL
            output$IDplot <- NULL
          }
        }


        output$selectSampleFromUploadedVCF <- renderUI(
          {

            sample.names <- colnames(counts.catalog[[1]])
            radioButtons(inputId = "sampleNameFromUploadedVCF",
                         label = "Select the sample from uploaded VCF",
                         choices = sample.names,
                         selected = character(0))
          }
        )

        output$spectraPlotFromVCF <- renderUI(
          {
            
            tabsetPanel(type = "tabs",
                        tabPanel("SBS96", plotOutput("SBS96plot")),
                        tabPanel("SBS192", plotOutput("SBS192plot")),
                        tabPanel("SBS1536", plotOutput("SBS1536plot")),
                        tabPanel("DBS78", plotOutput("DBS78plot")),
                        tabPanel("DBS136", plotOutput("DBS136plot")),
                        tabPanel("DBS144", plotOutput("DBS144plot")),
                        tabPanel("ID", plotOutput("IDplot"))
            )
          }
        )

        observeEvent(input$sampleNameFromUploadedVCF, {
          output$SBS96plot <- renderPlot({
            catSBS96 <-
              counts.catalog$catSBS96[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catSBS96)
          })


          output$SBS192plot <- renderPlot({
            catSBS192 <-
              counts.catalog$catSBS192[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catSBS192)
          })

          output$SBS1536plot <- renderPlot({
            catSBS1536 <-
              counts.catalog$catSBS1536[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catSBS1536)
          })

          output$DBS78plot <- renderPlot({
            catDBS78 <-
              counts.catalog$catDBS78[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catDBS78)
          })

          output$DBS136plot <- renderPlot({
            catDBS136 <-
              counts.catalog$catDBS136[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catDBS136)
          })

          output$DBS144plot <- renderPlot({
            catDBS144 <-
              counts.catalog$catDBS144[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catDBS144)
          })

        })

        output$selectSampleFromVCFForAttribution <- renderUI(
          {

            sample.names <- colnames(counts.catalog[[1]])
            selectInput(inputId = "selectedSampleFromVCFForAttribution",
                        label = "Select the sample from uploaded VCF",
                        choices = sample.names)
          }
        )

        observeEvent(input$selectedSampleFromVCFForAttribution, {
          output$selectCancerType <- renderUI(
            {
              cancer.types <-
                c(colnames(CancerTypeToExposureStatData()), "Unknown")
              selectInput(inputId = "selectedCancerType",
                          label = "Select the cancer type",
                          choices = cancer.types,
                          selected = "Biliary-AdenoCA")
            }
          )
        })

        observeEvent(input$selectedSampleFromVCFForAttribution, {
        output$choosecatalogtype <- renderUI(
          {
            catalog.type <- c("SBS96", "SBS192", "DBS78", "ID")
            selectInput(inputId = "selectedCatalogType",
                        label = "Select the catalog type",
                        choices = catalog.type,
                        selected = input.catalog.type)
          }
        )
        })

        observeEvent(input$selectedSampleFromVCFForAttribution, {
          output$chooseSigSubsetForSampleFromVCF <- renderUI(
            {
              sig1 <- PCAWG7::signature[["genome"]]
              sig2 <- sig1[[input$selectedCatalogType]]
              sig.universe <- colnames(sig2)


              foo <- CancerTypeToSigSubset(cancer.type = input$selectedCancerType,
                                           tumor.cohort = "PCAWG",
                                           sig.type = input$selectedCatalogType,
                                           region = "genome")
              selected.sig.universe <- colnames(foo)
              selectInput(inputId = "selectedSigSubset1",
                          label = paste0("These signatures were preselected based ",  
                                         "on cancer type. Move your cursor to click one ",
                                         "signature and press Backspace key to exclude ", 
                                         "the signature. Click the empty space inside ",
                                         "the box below to add new signature."),
                          choices = sig.universe,
                          selected = selected.sig.universe,
                          multiple = TRUE)
            }
          )
        }
        )

        observeEvent(input$selectedSampleFromVCFForAttribution, {
          output$analyzeButtonForVCF <- renderUI(
            {
              actionButton(inputId = "submitAttributionForVCF", label = "Analyze",
                           style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
            }
          )
        })
      } else if (input$vcftype == "strelka.id") {
        old.error.ids <- ids$error
        ids <<- ProcessStrelkaIDVCFs(input, output, file, ids)
        ids$error <- append(ids$error, old.error.ids)
      } else if (input$vcftype == "mutect") {
        old.error.ids <- ids$error
        ids <<- ProcessMutectVCFs(input, output, file, ids)
        ids$error <- append(ids$error, old.error.ids)
      }

      # When the downloadHandler function runs, increment rv$downloadFlag
      #rv$downloadFlag <- rv$downloadFlag + 1
    })
  
  ShowTwoButtons <- function() {
    output$showSpectra <- renderUI(
      actionButton(inputId = "showSpectraFromCatalog", label = "Show spectra",
                   style="color: #fff; background-color: #337ab7;
                          border-color: #2e6da4;")
    )
    output$sigAttribution <- renderUI(
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
      shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
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
            actionButton(inputId = "clickToSigAttribution",
                         label = "Signature attribution",
                         style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
          )
        }
      )
  })
  
  # When user clicks the action button on Show spectra page, direct user to the relevant tab
  observeEvent(input$clickToSigAttribution, {
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
    if (input.catalog.type == "SBS96") {
      output$SBS96plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])}, width = 800, height = 260)
    } else if (input.catalog.type == "SBS192") {
      output$SBS192plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input.catalog.type == "SBS1536") {
      output$SBS1536plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input.catalog.type == "DBS78") {
      output$DBS78plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input.catalog.type == "DBS136") {
      output$DBS136plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input.catalog.type == "DBS144") {
      output$DBS144plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input.catalog.type == "ID") {
      output$IDplot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    }
  })
  
  # When user selects the sample from uploaded catalog, show
  # the sample's mutational spectrum
  observeEvent(input$selectedSampleFromUploadedCatalog, {
    
    output$spectraPlotFromCatalog <- renderUI(
      {
        output$spectrum <- renderPlot(
          {
            PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                                drop = FALSE])
          }, width = 800, height = 260)
        plotOutput(outputId = "spectrum")
      }
    )
        
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
      
      output$selectSampleFromCatalogForAttribution <- renderUI(
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
          selectInput(inputId = "selectedSampleFromCatalogForAttribution",
                      label = "Select the sample from uploaded spectra",
                      choices = sample.names)
        }
      )
      
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
      
      output$selectCancerType <- renderUI(
        {
          cancer.types <-
            c("Unknown", colnames(CancerTypeToExposureStatData()))
          selectInput(inputId = "selectedCancerType",
                      label = "Select the cancer type",
                      choices = cancer.types)
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
      
      output$choosecatalogtype <- renderUI(
        {
          catalog.type <- c("SBS96", "SBS192", "DBS78", "ID")
          selectInput(inputId = "selectedCatalogType",
                      label = "Select the catalog type",
                      choices = catalog.type,
                      selected = input.catalog.type)
        }
      )
      
      output$uploadedCatalogType <- renderUI(
        {
            p(tags$b("Mutation type: "), input.catalog.type)
        }
      )
    })

  observeEvent(input$selectedCancerType2, {
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
    output$mytable <- DT::renderDataTable({
      DT::datatable(dat[sigsForAttribution(), ], escape = FALSE, rownames = FALSE) 
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
    
    # Hide the previous attribution results
    if (attribution.results == TRUE) {
      shinyjs::hide(id = "attributionResults")
    }
    
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
        
        if (retval$success == FALSE || is.null(retval$success)) {
          output$attributionResults <- renderUI({
            output$attributionMessage <- 
              renderText(paste0("The algorithm could not find the optimal number of ", 
                                "signatures that explain the spectrum. Please reduce the ", 
                                "number of signatures used."))
            tagList(
              textOutput(outputId = "attributionMessage")
            )
          })
        } else {
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
        }
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
    
    if (FALSE) {
      if (input.catalog.type == "SBS96") {
        output$sigTestButton2 <- renderUI(
          { #tags$head(tags$style(make_css(list('.btn', 'white-space', 'pre-wrap'))))
            actionButton(inputId = "submitSigTest2", 
                         label = "Check for artifacts/rare signatures",
                         style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
          })
      }
    }
    
    # Return something other than the future so we don't block the UI
    NULL
  })
  
  observeEvent(input$selectedSigSubset2, {
    output$analyzeButtonOnTop <- renderUI(
      {
        actionButton(inputId = "submitAttributionOnTop", label = "Analyze",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
      }
    )
    
    output$analyzeButton2 <- renderUI(
      {
        actionButton(inputId = "submitAttribution2", label = "Analyze",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
      }
    )
  })
  
  # Asynchronous programming code starts here
  submitAttribution <- reactive({
    list(input$submitAttributionOnTop, input$submitAttribution2)
  })
  
  observeEvent(
    submitAttribution(),
    {
      if(is.null(input$submitAttributionOnTop) && is.null(input$submitAttribution2)){
        return()
      }
      
      if(input$submitAttributionOnTop == 0 && input$submitAttribution2 == 0){
        return()
      }
      
      if(length(input$selectedSigSubset2) == 0) {
        showNotification(ui = "Error:", 
                         action = "No signatures selected for attribution analysis",
                         type = "error")
        return()
      }
      
      #Don't do anything if in the middle of a run
      if(running())
        return(NULL)
      running(TRUE)
      
      # Hide the previous attribution results
      if (attribution.results == TRUE) {
        shinyjs::hide(id = "attributionResults")
      }
      
      spect <- catalog[, input$selectedSampleFromCatalogForAttribution, drop = FALSE]
      catalog.type <- input$selectedCatalogType
      cancer.type <- input$selectedCancerType
      region <- input$region2
      
      if (catalog.type == "SBS192") {
        sig.universe <- 
          PCAWG7::signature[["genome"]][[catalog.type]][, input$selectedSigSubset2, drop = FALSE]
      } else if (catalog.type == "SBS96") {
        sig.universe <- 
          COSMIC.v3.genome.SBS96.sigs[, input$selectedSigSubset2, drop = FALSE]
      } else {
        sig.universe <- 
          PCAWG7::signature[[region]][[catalog.type]][, input$selectedSigSubset2, drop = FALSE]
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
          
          if (retval$success == FALSE || is.null(retval$success)) {
            output$attributionResults <- renderUI({
              output$attributionMessage <- 
                renderText(paste0("The algorithm could not find the optimal number of ", 
                                  "signatures that explain the spectrum. Please reduce the ", 
                                  "number of signatures used."))
              tagList(
                textOutput(outputId = "attributionMessage")
              )
            })
          } else {
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
            
            retval <- PrepareAttributionResults(input, output, input.catalog.type, 
                                                plotdata)
            attribution.results <<- retval$attribution.results
        }
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
      
      if (input.catalog.type == "SBS96") {
      output$sigTestButton2 <- renderUI(
        { #tags$head(tags$style(make_css(list('.btn', 'white-space', 'pre-wrap'))))
          actionButton(inputId = "submitSigTest2", 
                       label = "Check for artifacts/rare signatures",
                       style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4;padding:4px; ")
        })
      }
      
      # Return something other than the future so we don't block the UI
      NULL
    })
  
  # Synchronous programming code ends here
  
  # Send interrupt signal to future
  observeEvent(input$cancel,{
    if(running()) {
      interruptor$interrupt("Task cancelled")
    }
  })
  
  # Save extra values in state$values when we bookmark
  onBookmark(function(state) {
    state$values$previous.plotdata <- plotdata
    state$values$previous.input.catalog.type <- input.catalog.type
  })
  
  # Read values from state$values when we restore
  onRestore(function(state) {
    shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
    shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab2]')
    
    # Set the value of attribution.results to be TRUE, so when user submit
    # analysis again, hide the attribution results page
    attribution.results <<- TRUE
    
    PrepareAttributionResults(input = input, output = output, 
                              input.catalog.type = state$values$previous.input.catalog.type, 
                              plotdata = state$values$previous.plotdata)
  })
  
  observeEvent(input$submitAttributionForVCF, {
    output$cancelButton <- renderUI(
      actionButton(inputId = "cancel", label = "Cancel")
    )
  })
  
  observeEvent(input$submitAttribution2, {
    output$cancelButton <- renderUI(
      actionButton(inputId = "cancel", label = "Cancel")
    )
  })
  
  observeEvent(input$submitAttributionOnTop, {
    output$cancelButton <- renderUI(
      actionButton(inputId = "cancel", label = "Cancel")
    )
  })
  
  # Exclude the Analyze button from bookmarking
  
  setBookmarkExclude(c("submitAttributionForVCF", "submitAttribution2", 
                       "submitAttributionOnTop"))
  
  observeEvent(input$submitAttributionForVCF, {
    output$bookmarkButton <- renderUI(
      bookmarkButton(label = "Save the session",
                     title = "This application's state can be saved on the server for 7 days.")
    )
  })
  
  observeEvent(input$submitAttribution2, {
    output$bookmarkButton <- renderUI(
      bookmarkButton(label = "Save the session",
                     title = "This application's state can be saved on the server for 7 days.")
    )
  })
  
  observeEvent(input$submitAttributionOnTop, {
    output$bookmarkButton <- renderUI(
      bookmarkButton(label = "Save the session",
                     title = "This application's state can be saved on the server for 7 days.")
    )
  })
  
  observeEvent(input$submitAttributionOnTop, {
    output$downloadResults <- renderUI({
      downloadButton(outputId = 'downloadAttributionResults', 
                     label = 'Download attribution results')
    })
  })
  
  observeEvent(input$submitAttribution2, {
    output$downloadResults <- renderUI({
      downloadButton(outputId = 'downloadAttributionResults', 
                     label = 'Download attribution results')
    })
  })
  
  # Show modal when button is clicked.
  observeEvent(input$submitSigTest2, {
    showModal(ui =  modalDialog(
      selectInput(inputId = "selectedArtifact",
                  label = "Select artifact signature to test",
                  choices = c("None", mSigAct::PossibleArtifacts())),
      
      selectInput(inputId = "selectedRareSig",
                  label = "Select rare signature to test",
                  choices = c("None", mSigAct::RareSignatures())),
      
      
      footer = tagList(
        modalButton(label = "Close"),
        actionButton(inputId = "testSigPresence", label = "Test",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
      )))
    
  })
  
  observeEvent(input$testSigPresence, {
    #hideTab(inputId = "attributionTabSet", target = "sigPresenceTestTab")
    removeTab(inputId = "attributionTabSet", target = "sigPresenceTestTab")
    withProgress(message = "Testing in progress", value = 0, detail="0%", {
      
      if (input$selectedArtifact == "None" && input$selectedRareSig == "None") {
        showNotification(ui = "Error:", 
                         action = "No signature selected for testing",
                         type = "error")
        return()
      } 
      
      sigs <- plotdata$sig.universe
      map.sigs.names <- plotdata$QP.best.MAP.exp$sig.id
      map.sigs <- sigs[, map.sigs.names, drop = FALSE]
      
      ConvertTextToLinks <- function(text, urls) {
        urls1 <- urls[text, ]
        refs <- paste0("<a href='",  urls1, "' target='_blank'>", 
                       text, "</a>")
        return(refs)
      }
      
      map.sig.refs <- ConvertTextToLinks(text = map.sigs.names, 
                                         urls = COSMIC.v3.SBS.sig.links)
      
      if (input$selectedArtifact != "None") {
        
        shinyjs::hide(id = "artifactSigTest")
        shinyjs::hide(id = "rareSigTest")
        artifact.sig.to.test <- 
          COSMIC.v3.genome.SBS96.sigs[, input$selectedArtifact, drop = FALSE]
        artifact.sigs <- cbind(artifact.sig.to.test, map.sigs)
        incProgress(0.1, detail = "Testing the artifact signature")
        artifact.sig.test <- 
          mSigAct::SignaturePresenceTest1(spectrum = plotdata$spect,
                                          sigs = artifact.sigs,
                                          target.sig.index = 1,
                                          m.opts = mSigAct::DefaultManyOpts())
        artifact.sig.refs.model1 <- ConvertTextToLinks(text = colnames(artifact.sigs), 
                                                       urls = COSMIC.v3.SBS.sig.links)
        
        # Change the first URL link to red color
        artifact.sig.refs.model1[1] <- 
          gsub(pattern = "target='_blank'", 
               replacement = "target='_blank' style = 'color: red'",
               x = artifact.sig.refs.model1[1])
        
        artifact.sig.test.model1 <- 
          data.frame(artifact.signature.presence.test = 
                       paste(artifact.sig.refs.model1, collapse = ","), 
                     log.likelihood = artifact.sig.test$with,
                     test.statistic = artifact.sig.test$statistic,
                     p.value = artifact.sig.test$chisq.p)
        artifact.sig.test.model2 <-
          data.frame(artifact.signature.presence.test = 
                       paste(map.sig.refs, collapse = ","), 
                     log.likelihood = artifact.sig.test$without)
        
        artifact.sig.test.output <- dplyr::bind_rows(artifact.sig.test.model1,
                                                     artifact.sig.test.model2)
        incProgress(0.4, detail = "Testing the rare signature")
      }
      
      if (input$selectedRareSig != "None") {
        shinyjs::hide(id = "artifactSigTest")
        shinyjs::hide(id = "rareSigTest")
        rare.sig.to.test <- 
          COSMIC.v3.genome.SBS96.sigs[, input$selectedRareSig, drop = FALSE]
        rare.sigs <- cbind(rare.sig.to.test, map.sigs)
        rare.sig.test <-
          mSigAct::SignaturePresenceTest1(spectrum = plotdata$spect,
                                          sigs = rare.sigs,
                                          target.sig.index = 1,
                                          m.opts = mSigAct::DefaultManyOpts())
        
        incProgress(0.4, detail = "Preparing output")
        rare.sig.refs.model1 <- ConvertTextToLinks(text = colnames(rare.sigs), 
                                                   urls = COSMIC.v3.SBS.sig.links)
        rare.sig.refs.model1[1] <- 
          gsub(pattern = "target='_blank'", 
               replacement = "target='_blank' style = 'color: red'",
               x = rare.sig.refs.model1[1])
        
        rare.sig.test.model1 <- 
          data.frame(rare.signature.presence.test = 
                       paste(rare.sig.refs.model1, collapse = ","), 
                     log.likelihood = rare.sig.test$with,
                     test.statistic = rare.sig.test$statistic,
                     p.value = rare.sig.test$chisq.p)
        rare.sig.test.model2 <-
          data.frame(rare.signature.presence.test = 
                       paste(map.sig.refs, collapse = ","), 
                     log.likelihood = rare.sig.test$without)
        
        rare.sig.test.output <- dplyr::bind_rows(rare.sig.test.model1,
                                                 rare.sig.test.model2)
      }
      
      #incProgress(0.4, detail = "Preparing output")
      
      if (exists("artifact.sig.test.output")) {
        output$artifactSigTest <- renderTable({
          artifact.sig.test.output
        }, sanitize.text.function = function(x) x, digits = 5)
      } else {
        output$artifactSigTest <- NULL
      }
      
      if (exists("rare.sig.test.output")) {
        output$rareSigTest <- renderTable({
          rare.sig.test.output
        }, sanitize.text.function = function(x) x, digits = 5)
      } else {
        output$rareSigTest <- NULL
      }

      
      appendTab(inputId = "attributionTabSet", 
                tab = tabPanel(title = "Signature presence test", 
                               uiOutput(outputId = "artifactSigTest"),
                               br(),
                               uiOutput(outputId = "rareSigTest"),
                               value = "sigPresenceTestTab"), 
                select = TRUE)
      
      setProgress(value = 1, detail = "Output is ready")
      #incProgress(0.1, detail = "Output is ready")
    })
    
    
    removeModal()
    
  })
  
  # When user clicks the "Remove notifications" button, all the previous
  # notifications(error, warning or message) will be removed
  observeEvent(input$remove, {
    RemoveAllNotifications(ids)
  })
  
  observeEvent(input$remove2, {
    RemoveAllNotifications(ids)
  })

}
