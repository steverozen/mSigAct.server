#' @import shiny
#' @import shinydashboard
app_server <- function(input, output,session) {
  # List the first level callModules here
  
  plotdata <- reactiveValues(spect = NULL, reconstructed.catalog = NULL,
                             sig.universe = NULL, QP.best.MAP.exp = NULL)
  
  hideTab(inputId = "panels", target = "tab4")
  
  hideTab(inputId = "panels", target = "tab5")
  
  # When user clicks the action link on Home page, direct user to the relevant tab
  observeEvent(input$linkToTab2, {
    shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                   selected = "tab2")
  })
  
  observeEvent(input$linkToTab3, {
    shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                   selected = "tab3")
  })

  # Create an empty list which can be used to store notification ids later
  ids <- list("error" = character(0), "warning" = character(0),
             "message" = character(0))

  # Create an empty list which can be used to store the return value of processing VCFs
  retval <- list()

  # Create a variable which can be used to store the uploaded catalog later
  catalog <- NA
  
  input.catalog.type <- NA
  
  plot.names <- vector(mode = "character")

  # Create reactiveValues object
  # and set flag to 0 to prevent errors with adding NULL
  #rv <- reactiveValues(downloadFlag = 0)

  # Download sample VCFs when user clicks the button
  output$downloadsampleVCFs <- downloadHandler(
    filename = function() {
      "sample-VCFs.zip"
    },
    content = function(file) {
      PrepareSampleVCFs(file)
    })

  # Download sample catalogs when user clicks the button
  output$downloadSampleSpectra <- downloadHandler(
    filename = function() {
      "mSigAct-sample-spectra.zip"
    },
    content = function(file) {
      PrepareSampleSpectra(file)
    })

  # Run ICAMS on sample Strelka SBS VCFs when user clicks the button
  output$runstrelkasbsvcfs <- downloadHandler(
    filename = function() {
      "ICAMS-test-run-Strelka-SBS-VCFs.zip"
    },
    content = function(file) {
      ids <<-
        RunICAMSOnSampleStrelkaSBSVCFs(output, file, ids)
    })

  # Run ICAMS on sample Mutect VCFs when user clicks the button
  output$runmutectvcfs <- downloadHandler(
    filename = function() {
      "ICAMS-test-run-Mutect-VCFs.zip"
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
  
  
  vcftype <- reactive({
    validate(
      need(input$vcftype %in% c("strelka.sbs", "strelka.id", "mutect"),
           label = "variant caller")
    )
    paste0(input$vcftype, "test")
  })
  
  output$testoutput <- renderText(
    vcftype()
  )
  
  # When user submits VCF to analysis, then the program will try to
  # generate a zip archive based on the input files and parameters
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$zipfile.name, ".zip")
    },

    content = function(file) {
      if (input$vcftype == "strelka.sbs") {
        # Generate a zip archive from Strelka SBS VCFs and
        # update the notification ids for errors, warnings
        # and messages

        result <- ProcessStrelkaSBSVCFs(input, output, file, ids)
        retval <<- result$retval
        ids <<- result$ids

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
          output$analyzeButton1 <- renderUI(
            {
              actionButton(inputId = "submitAttribution1", label = "Analyze",
                           style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
            }
          )
        })
      } else if (input$vcftype == "strelka.id") {
        ids <<- ProcessStrelkaIDVCFs(input, output, file, ids)
      } else if (input$vcftype == "mutect") {
        ids <<- ProcessMutectVCFs(input, output, file, ids)
      }

      # When the downloadHandler function runs, increment rv$downloadFlag
      #rv$downloadFlag <- rv$downloadFlag + 1
    })

  # When user submit uploaded spectra for analysis, create radio buttons for
  # user to select the sample
  observeEvent(input$submitSpectra, {
    showTab(inputId = "panels", target = "tab4")
    showTab(inputId = "panels", target = "tab5")
    output$selectSampleFromUploadedCatalog <-
      renderUI(
        {
          # catalog.info is a data frame that contains one row for each uploaded file,
          # and four columns "name", "size", "type" and "datapath".
          # "name": The filename provided by the web browser.
          # "size": The size of the uploaded data, in bytes.
          # "type": The MIME type reported by the browser.
          # "datapath": The path to a temp file that contains the data that was uploaded.
          catalog.info <- input$upload.spectra
          catalog.paths <- catalog.info$datapath
          uploaded.catalog <- ICAMS::ReadCatalog(file = catalog.paths,
                                                 ref.genome = input$ref.genome2,
                                                 region = input$region2)
          catalog <<- uploaded.catalog
          
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
          radioButtons(inputId = "selectedSampleFromUploadedCatalog",
                       label = "Select the sample from uploaded catalog",
                       choices = sample.names,
                       selected = character(0))
        }
      )
  })

  # When user submit new catalog for analysis, remove the previous plots
  observeEvent(input$submitSpectra, {
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

  observeEvent(input$submitSpectra, {
    output$selectSampleFromCatalogForAttribution <- renderUI(
      { 
        # catalog.info is a data frame that contains one row for each uploaded file,
        # and four columns "name", "size", "type" and "datapath".
        # "name": The filename provided by the web browser.
        # "size": The size of the uploaded data, in bytes.
        # "type": The MIME type reported by the browser.
        # "datapath": The path to a temp file that contains the data that was uploaded.
        catalog.info <- input$upload.spectra
        catalog.paths <- catalog.info$datapath
        uploaded.catalog <- ICAMS::ReadCatalog(file = catalog.paths,
                                               ref.genome = input$ref.genome2,
                                               region = input$region2)
        catalog <<- uploaded.catalog
        
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
  })

  observeEvent(input$submitSpectra, {
  output$selectCancerType <- renderUI(
    {
      cancer.types <-
        c("Unknown", colnames(CancerTypeToExposureStatData()))
      selectInput(inputId = "selectedCancerType",
                  label = "Select the cancer type",
                  choices = cancer.types)
    }
  )
  })

  observeEvent(input$submitSpectra, {
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
            selected.sig.universe <- colnames(tmp)
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
      output$analyzeButton2 <- renderUI(
        {
          actionButton(inputId = "submitAttribution2", label = "Analyze",
                       style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
        }
      )
  })
  
  observeEvent(input$selectedSigSubset2, {
    output$analyzeButton3 <- renderUI(
      {
        actionButton(inputId = "submitAttribution3", label = "Analyze",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
      }
    )
  })
  
  observeEvent(input$submitAttribution2, {
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
      sig.universe <- PCAWG7::signature[["genome"]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- PCAWG7::exposure.stats$PCAWG[["SBS96"]][[cancer.type]][colnames(sig.universe), ]
    } else {
      sig.universe <- PCAWG7::signature[[region]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- PCAWG7::exposure.stats$PCAWG[[catalog.type]][[cancer.type]][colnames(sig.universe), ]
    }
    
    sig.names <- rownames(sigs.prop)
    sigs.prop <- unlist(sigs.prop[ , 2])
    names(sigs.prop) <- sig.names
    
    if (FALSE) {
      QP.exposure <- GetExposureWithConfidence(catalog = spect,
                                               sig.universe = sig.universe,
                                               num.of.bootstrap.replicates = 10000,
                                               method = decomposeQP,
                                               conf.int = 0.95)
      
      updated.sig.universe <- sig.universe[ , rownames(QP.exposure)]
      updated.sigs.prop <- sigs.prop[rownames(QP.exposure)]
    }
    
    mapout <-
      mSigAct::MAPAssignActivity1(
        spect = spect,
        sigs = sig.universe,
        sigs.presence.prop = sigs.prop,
        max.level = length(sigs.prop) - 1,
        p.thresh = 0.01,
        eval_f = mSigAct::ObjFnBinomMaxLHNoRoundOK,
        m.opts = mSigAct::DefaultManyOpts(),
        max.mc.cores = 100)
    
    MAP.best.exp <- mapout$MAP
    
    QP.exp <- 
      mSigAct:::OptimizeExposureQP(spect, sig.universe[ , 
                                                        MAP.best.exp$sig.id, 
                                                        drop = FALSE])
    QP.best.MAP.exp <-
      tibble::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)
    
    r.qp <- mSigAct::ReconstructSpectrum(sig.universe, exp = QP.exp, use.sig.names = TRUE)
    reconstructed.catalog0 <- as.catalog(r.qp, ref.genome = input$ref.genome2, 
                                         region = input$region2)
    
    cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
    
    colnames(reconstructed.catalog0) <- 
      paste0("reconstructed (cosine similarity = ", cossim, ")")
    reconstructed.catalog <- round(reconstructed.catalog0)
    
    
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
    
    for (i in length(plot.names)) {
      shinyjs::show(id = plot.names[i])
    }
    
    
  })
  
  observeEvent(input$submitAttribution3, {
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
      sig.universe <- PCAWG7::signature[["genome"]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- PCAWG7::exposure.stats$PCAWG[["SBS96"]][[cancer.type]][colnames(sig.universe), ]
    } else {
      sig.universe <- PCAWG7::signature[[region]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- PCAWG7::exposure.stats$PCAWG[[catalog.type]][[cancer.type]][colnames(sig.universe), ]
    }
    
    sig.names <- rownames(sigs.prop)
    sigs.prop <- unlist(sigs.prop[ , 2])
    names(sigs.prop) <- sig.names
    
    if (FALSE) {
      QP.exposure <- GetExposureWithConfidence(catalog = spect,
                                               sig.universe = sig.universe,
                                               num.of.bootstrap.replicates = 10000,
                                               method = decomposeQP,
                                               conf.int = 0.95)
      
      updated.sig.universe <- sig.universe[ , rownames(QP.exposure)]
      updated.sigs.prop <- sigs.prop[rownames(QP.exposure)]
    }
    
    mapout <-
      mSigAct::MAPAssignActivity1(
        spect = spect,
        sigs = sig.universe,
        sigs.presence.prop = sigs.prop,
        max.level = length(sigs.prop) - 1,
        p.thresh = 0.01,
        eval_f = mSigAct::ObjFnBinomMaxLHNoRoundOK,
        m.opts = mSigAct::DefaultManyOpts(),
        max.mc.cores = 100)
    
    MAP.best.exp <- mapout$MAP
    
    QP.exp <- 
      mSigAct:::OptimizeExposureQP(spect, sig.universe[ , 
                                                        MAP.best.exp$sig.id, 
                                                        drop = FALSE])
    QP.best.MAP.exp <-
      tibble::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)
    
    r.qp <- mSigAct::ReconstructSpectrum(sig.universe, exp = QP.exp, use.sig.names = TRUE)
    reconstructed.catalog0 <- as.catalog(r.qp, ref.genome = input$ref.genome2, 
                                         region = input$region2)
    
    cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
    
    colnames(reconstructed.catalog0) <- 
      paste0("reconstructed (cosine similarity = ", cossim, ")")
    reconstructed.catalog <- round(reconstructed.catalog0)
    
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
    
    for (i in length(plot.names)) {
      shinyjs::show(id = plot.names[i])
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
  
  # Exclude the Analyze button from bookmarking
  
  setBookmarkExclude(c("submitAttribution1", "submitAttribution2", 
                       "submitAttribution3"))
  
  observeEvent(input$submitAttribution1, {
    output$bookmarkButton <- renderUI(
      bookmarkButton()
    )
  })
  
  observeEvent(input$submitAttribution2, {
    output$bookmarkButton <- renderUI(
      bookmarkButton()
    )
  })
  
  observeEvent(input$submitAttribution3, {
    output$bookmarkButton <- renderUI(
      bookmarkButton()
    )
  })
  
  # When user clicks the "Remove notifications" button, all the previous
  # notifications(error, warning or message) will be removed
  observeEvent(input$remove, {
    RemoveAllNotifications(ids)
  })

}
