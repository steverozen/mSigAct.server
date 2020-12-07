# Cannot use plan(multicore), otherwise the progress bar for asynchronous
# process will not work properly
future::plan(future::multisession)

#' @import tibble
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
  
  output$exampleSpectraText <- renderText("Load example spectra")
  
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
  
  catalog.path <- input.catalog.type <- NA
  
  showSBS192Catalog <- TRUE
  
  plot.names <- vector(mode = "character")

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
                          border-color: #2e6da4")
    )
    output$sigAttribution <- renderUI(
      actionButton(inputId = "sigAttributionFromCatalog", 
                   label = "Signature attribution",
                   style="color: #fff; background-color: #337ab7;
                          border-color: #2e6da4"))
  }
  
  observeEvent(input$preloadSBS96Spectra, {
    
    shinyWidgets::updatePickerInput(session = session,
                                    inputId = "ref.genome2",
                                    selected = "hg19")
    shinyWidgets::updatePickerInput(session = session,
                                    inputId = "region2",
                                    selected = "genome")
    catalog.path <<- system.file("extdata/SBS96-mSigAct-example-spectra.csv", 
                                 package = "ICAMS.shiny")
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
                                 package = "ICAMS.shiny")
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
                                 package = "ICAMS.shiny")
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
                                 package = "ICAMS.shiny")
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
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
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
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab")
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
                                   selected = "sigAttributionTab")
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
                      label = "Select the sample from uploaded catalog",
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


  observeEvent(input$selectedCancerType, {
    if (input$selectedCancerType != "Unknown") {
      output$chooseSigSubsetForSampleFromCatalog <- renderUI(
        {
          sig1 <- PCAWG7::signature[["genome"]]
          sig2 <- sig1[[input$selectedCatalogType]]
          sig.universe <- colnames(sig2)
          
          if (input$selectedCancerType == "Unknown") {
            selected.sig.universe <- NULL
          } else {
            tmp <- CancerTypeToSigSubset(cancer.type = input$selectedCancerType,
                                         tumor.cohort = "PCAWG",
                                         sig.type = input$selectedCatalogType,
                                         region = "genome")
            selected.sig.universe0 <- colnames(tmp)
            
            # Exclude possible artifact signatures
            # possible.artifacts <- mSigAct::PossibleArtifacts()
            possible.artifacts <- NULL
            
            selected.sig.universe <- 
              setdiff(selected.sig.universe0, possible.artifacts)
          }
          
          selectInput(inputId = "selectedSigSubset2",
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
    
  })

  observeEvent(input$selectedSigSubset2, {
    output$analyzeButtonOnTop <- renderUI(
      {
        actionButton(inputId = "submitAttributionOnTop", label = "Analyze",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
      }
    )
    
    output$analyzeButton2 <- renderUI(
      {
        actionButton(inputId = "submitAttribution2", label = "Analyze",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
      }
    )
  })
  
  # Synchronous programming code starts here
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
      
      #Don't do anything if in the middle of a run
      if(running())
        return(NULL)
      running(TRUE)
      
      
      if (length(plot.names) > 0) {
        for (i in 1:length(plot.names)) {
          shinyjs::hide(id = plot.names[i])
        }
      }
      
      spect <- catalog[, input$selectedSampleFromCatalogForAttribution, drop = FALSE]
      catalog.type <- input$selectedCatalogType
      cancer.type <- input$selectedCancerType
      region <- input$region2
      
      if (catalog.type == "SBS192") {
        sig.universe <- 
          PCAWG7::signature[["genome"]][[catalog.type]][, input$selectedSigSubset2]
        sigs.prop <- 
          PCAWG7::exposure.stats$PCAWG[["SBS96"]][[cancer.type]][colnames(sig.universe), ]
      } else {
        sig.universe <- 
          PCAWG7::signature[[region]][[catalog.type]][, input$selectedSigSubset2]
        sigs.prop <- 
          PCAWG7::exposure.stats$PCAWG[[catalog.type]][[cancer.type]][colnames(sig.universe), ]
      }
      
      sig.names <- rownames(sigs.prop)
      sigs.prop <- unlist(sigs.prop[ , 2])
      names(sigs.prop) <- sig.names
      
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
            
            # TODO: Need to change the callback function in mSigAct::MAPAssignActivityInternal
            # each.level.callback.fn(
            # value = 1/max.level - 0.01,
            # detail = paste0("Testing removal of subsets of ", df, " signatures (",
            #                 length(subsets2), " subsets)"))
            progress$inc(amount = value, detail = detail)
            interruptor$execInterrupts()
          }
          
          retval <- mSigAct::MAPAssignActivity1(
            spect = spect,
            sigs = sig.universe,
            sigs.presence.prop = sigs.prop,
            max.level = length(sigs.prop) - 1,
            p.thresh = 0.01,
            eval_f = mSigAct::ObjFnBinomMaxLHRound,
            eval_g_ineq = mSigAct::g_ineq_for_ObjFnBinomMaxLH2,
            m.opts = mSigAct::DefaultManyOpts(),
            max.mc.cores = 50,
            progress.monitor = updateProgress
          )
          
          return(retval)
        }, seed = TRUE) %...>% {
          retval <- .
          
          
          if (retval$success == FALSE || is.null(retval$success)) {
            output$sigContributionPlot <- renderUI({
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
              tibble::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)
            
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
            
            max_plots <- nrow(QP.best.MAP.exp) + 2
            output$sigContributionPlot <- renderUI({
              plot_output_list <- lapply(1:max_plots, function(i) {
                plotname <- paste("plot", i, sep="")
                plot.names[i] <<- plotname
                plotOutput(plotname)
              })
              
              tagList(plot_output_list)
            })
            
            for (i in 1:max_plots) {
              # Need local so that each item gets its own number. Without it, the value
              # of i in the renderPlot() will be the same across all instances, because
              # of when the expression is evaluated.
              local({
                my_i <- i
                plotname <- paste("plot", my_i, sep="")
                
                if (my_i == 1) {
                  output[[plotname]] <- renderPlot(
                    expr = ICAMS::PlotCatalog(spect)
                    , width = 800, height = 260
                  )
                } else if (my_i == 2) {
                  output[[plotname]] <- renderPlot(
                    expr = ICAMS::PlotCatalog(reconstructed.catalog)
                    , width = 800, height = 260
                  )
                } else {
                  output[[plotname]] <- renderPlot({
                    sig.name <- QP.best.MAP.exp$sig.id[my_i-2]
                    sig.catalog <- sig.universe[, sig.name, drop = FALSE]
                    colnames(sig.catalog) <- 
                      paste0(sig.name, " (exposure = ", 
                             round(QP.best.MAP.exp$QP.best.MAP.exp[my_i-2]), ")")
                    ICAMS::PlotCatalog(sig.catalog)
                    
                  }, width = 800, height = 260)
                }
              })
            }
            
            for (i in length(plot.names)) {
              shinyjs::show(id = plot.names[i])
            }
            
          }
          
        } %...>% {
          PrepareAttributionResults(input, output, input.catalog.type, 
                                    file, plotdata)
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
  
  # Synchronous programming code ends here
  
  # Send interrupt signal to future
  observeEvent(input$cancel,{
    if(running()) {
      interruptor$interrupt("Task cancelled")
    }
  })
  
  # Save extra values in state$values when we bookmark
  onBookmark(function(state) {
    state$values$previousSpect <- plotdata$spect
    state$values$previousReconstructed.catalog <- plotdata$reconstructed.catalog
    state$values$previousSig.universe <- plotdata$sig.universe
    state$values$previousQP.best.MAP.exp <- plotdata$QP.best.MAP.exp
  })
  
  # Read values from state$values when we restore
  onRestore(function(state) {
    spect <- state$values$previousSpect
    reconstructed.catalog <- state$values$previousReconstructed.catalog
    sig.universe <- state$values$previousSig.universe
    QP.best.MAP.exp <- state$values$previousQP.best.MAP.exp
    
    #spect <- plotdata$spect
    #reconstructed.catalog <- plotdata$reconstructed.catalog 
    #sig.universe <- plotdata$sig.universe 
    #QP.best.MAP.exp <- plotdata$QP.best.MAP.exp 
    
    max_plots <- nrow(QP.best.MAP.exp) + 2
    output$sigContributionPlot <- renderUI({
      plot_output_list <- lapply(1:max_plots, function(i) {
        plotname <- paste("plot", i, sep="")
        plot.names[i] <<- plotname
        plotOutput(plotname)
      })
      
      tagList(plot_output_list)
    })
    
    for (i in 1:max_plots) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_i <- i
        plotname <- paste("plot", my_i, sep="")
        
        if (my_i == 1) {
          output[[plotname]] <- renderPlot(
            expr = ICAMS::PlotCatalog(spect)
            #width = 800, height = 200, 
          )
        } else if (my_i == 2) {
          output[[plotname]] <- renderPlot(
            expr = ICAMS::PlotCatalog(reconstructed.catalog)
            #width = 800, height = 200)
          )
        } else {
          output[[plotname]] <- renderPlot({
            sig.name <- QP.best.MAP.exp$sig.id[my_i-2]
            sig.catalog <- sig.universe[, sig.name, drop = FALSE]
            colnames(sig.catalog) <- 
              paste0(sig.name, " (exposure = ", 
                     round(QP.best.MAP.exp$QP.best.MAP.exp[my_i-2]), ")")
            ICAMS::PlotCatalog(sig.catalog)
            
          }) #width = 800, height = 200)
        }
      })
    }
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
                     title = "This application's state can be saved on the server for 12 hours.")
    )
  })
  
  observeEvent(input$submitAttribution2, {
    output$bookmarkButton <- renderUI(
      bookmarkButton(label = "Save the session",
                     title = "This application's state can be saved on the server for 12 hours.")
    )
  })
  
  observeEvent(input$submitAttributionOnTop, {
    output$bookmarkButton <- renderUI(
      bookmarkButton(label = "Save the session",
                     title = "This application's state can be saved on the server for 12 hours.")
    )
  })
  
  observeEvent(input$submitAttributionOnTop, {
    output$attributionResults <- renderUI({
      downloadButton(outputId = 'downloadAttributionResults', 
                     label = 'Download attribution results')
    })
  })
  
  observeEvent(input$submitAttribution2, {
    output$attributionResults <- renderUI({
      downloadButton(outputId = 'downloadAttributionResults', 
                     label = 'Download attribution results')
    })
  })
  
  # Download attribution results when user clicks the button
  output$downloadAttributionResults <- downloadHandler(
    filename = function() {
      "mSigAct-signature-attribution-results.pdf"
    },
    content = function(file) {
      PrepareAttributionResults(input, file, plotdata)
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
