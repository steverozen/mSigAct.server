#' @import shiny
app_server <- function(input, output,session) {
  # List the first level callModules here

  # Create an empty list which can be used to store notification ids later
  ids <- list("error" = character(0), "warning" = character(0),
             "message" = character(0))

  # Create an empty list which can be used to store the return value of processing VCFs
  retval <- list()

  # Create a variable which can be used to store the uploaded catalog later
  catalog <- NA

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
  output$downloadSampleCatalogs <- downloadHandler(
    filename = function() {
      "sample-catalogs.zip"
    },
    content = function(file) {
      PrepareSampleCatalogs(file)
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
          output$selectcancertype <- renderUI(
            {
              cancer.types <-
                c(colnames(CancerTypeToExposureStatData()), "Unknown")
              selectInput(inputId = "selectedcancertype",
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
                        selected = "SBS96")
          }
        )
        })

        observeEvent(input$selectedSampleFromVCFForAttribution, {
          output$chooseSigSubsetForSampleFromVCF <- renderUI(
            {
              sig1 <- PCAWG7::signature[["genome"]]
              sig2 <- sig1[[input$selectedCatalogType]]
              sig.universe <- colnames(sig2)


              foo <- CancerTypeToSigSubset(cancer.type = input$selectedcancertype,
                                           tumor.cohort = "PCAWG",
                                           sig.type = input$selectedCatalogType,
                                           region = "genome")
              selected.sig.universe <- colnames(foo)
              selectInput(inputId = "selectedSigSubset1",
                          label = paste0("The following signatures are preselected according ",  
                                         "to previous assignment. Click the empty space ",
                                         "inside the box below to add new signature or use ",
                                         "Backspace to exclude signature"),
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

  # When user submit uploaded catalog for analysis, create radio buttons for
  # user to select the sample
  observeEvent(input$submitCatalog, {
    output$selectSampleFromUploadedCatalog <-
      renderUI(
        {
          # catalog.info is a data frame that contains one row for each uploaded file,
          # and four columns "name", "size", "type" and "datapath".
          # "name": The filename provided by the web browser.
          # "size": The size of the uploaded data, in bytes.
          # "type": The MIME type reported by the browser.
          # "datapath": The path to a temp file that contains the data that was uploaded.
          catalog.info <- input$upload.catalogs
          catalog.paths <- catalog.info$datapath
          uploaded.catalog <- ICAMS::ReadCatalog(file = catalog.paths,
                                                 ref.genome = input$ref.genome2,
                                                 region = input$region2)
          catalog <<- uploaded.catalog

          sample.names <- colnames(catalog)
          radioButtons(inputId = "selectedSampleFromUploadedCatalog",
                       label = "Select the sample from uploaded catalog",
                       choices = sample.names,
                       selected = character(0))
        }
      )
  })

  # When user submit new catalog for analysis, remove the previous plots
  observeEvent(input$submitCatalog, {
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

    if (input$catalogType == "SBS96") {
      output$SBS96plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input$catalogType == "SBS192") {
      output$SBS192plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input$catalogType == "SBS1536") {
      output$SBS1536plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input$catalogType == "DBS78") {
      output$DBS78plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input$catalogType == "DBS136") {
      output$DBS136plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input$catalogType == "DBS144") {
      output$DBS144plot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    } else if (input$catalogType == "ID") {
      output$IDplot <- renderPlot({
        PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog,
                            drop = FALSE])})
    }
  })

  observeEvent(input$submitCatalog, {
    output$selectSampleFromCatalogForAttribution <- renderUI(
      {
        sample.names <- colnames(catalog)
        selectInput(inputId = "selectedSampleFromCatalogForAttribution",
                    label = "Select the sample from uploaded catalog",
                    choices = sample.names)
      }
    )
  })

  observeEvent(input$submitCatalog, {
  output$selectcancertype <- renderUI(
    {
      cancer.types <-
        c(colnames(CancerTypeToExposureStatData()), "Unknown")
      selectInput(inputId = "selectedcancertype",
                  label = "Select the cancer type",
                  choices = cancer.types,
                  selected = "Biliary-AdenoCA")
    }
  )
  })

  observeEvent(input$submitCatalog, {
  output$choosecatalogtype <- renderUI(
    {
      catalog.type <- c("SBS96", "SBS192", "DBS78", "ID")
      selectInput(inputId = "selectedCatalogType",
                  label = "Select the catalog type",
                  choices = catalog.type,
                  selected = "SBS96")
    }
  )
  })

  observeEvent(input$selectedSampleFromCatalogForAttribution, {
    output$chooseSigSubsetForSampleFromCatalog <- renderUI(
      {
        sig1 <- PCAWG7::signature[["genome"]]
        sig2 <- sig1[[input$selectedCatalogType]]
        sig.universe <- colnames(sig2)


        foo <- CancerTypeToSigSubset(cancer.type = input$selectedcancertype,
                                     tumor.cohort = "PCAWG",
                                     sig.type = input$selectedCatalogType,
                                     region = "genome")
        selected.sig.universe <- colnames(foo)
        selectInput(inputId = "selectedSigSubset2",
                    label = paste0("The following signatures are preselected according ",  
                                   "to previous assignment. Click the empty space ",
                                   "inside the box below to add new signature or use ",
                                   "Backspace to exclude signature"),
                    choices = sig.universe,
                    selected = selected.sig.universe,
                    multiple = TRUE)
      }
    )
  })

  observeEvent(input$selectedSampleFromCatalogForAttribution, {
    output$analyzeButton2 <- renderUI(
      {
        actionButton(inputId = "submitAttribution2", label = "Analyze",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4")
      }
    )
  })

  observeEvent(input$submitAttribution2, {
    
    spect <- catalog[, input$selectedSampleFromCatalogForAttribution, drop = FALSE]
    catalog.type <- input$selectedCatalogType
    cancer.type <- input$selectedcancertype
    region <- input$region2
    
    if (catalog.type == "SBS192") {
      sig.universe <- PCAWG7::signature[["genome"]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- PCAWG7::exposure.stats$PCAWG[["SBS96"]][[cancer.type]]
    } else {
      sig.universe <- PCAWG7::signature[[region]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- PCAWG7::exposure.stats$PCAWG[[catalog.type]][[cancer.type]]
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
    
    xx <- mSigAct:::ListOfList2Tibble(mapout)
    
    best <- dplyr::arrange(xx, .data$MAP)[nrow(xx),  ]
    names.best <- names(best[["exp"]])
    best.exp <- best[["exp"]][[1]]
    if (is.null(names(best.exp))) {
      names(best.exp) <- names.best
    }
    MAP.best.exp <- tibble::tibble(sig.id = names(best.exp), best.exp )
    
    QP.exp <- mSigAct:::OptimizeExposureQP(spect, sig.universe[ , names(best.exp), drop = FALSE])
    QP.best.MAP.exp <-
      tibble::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)
    
    r.qp <- mSigAct::ReconstructSpectrum(sig.universe, exp = QP.exp, use.sig.names = TRUE)
    reconstructed.catalog <- as.catalog(r.qp)
    
    cossim <- round(mSigAct::cossim(spect, reconstructed.catalog), 5)
    
    colnames(reconstructed.catalog) <- paste0("MAP+QP ", cossim)
    reconstructed.catalog1 <- round(reconstructed.catalog)
    
    
    if (input$selectedCatalogType == "SBS96") {
      output$SBS96SpectrumPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(spect),
        width = 800, height = 200
      )
      output$SBS96AttributionPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(reconstructed.catalog1),
        width = 800, height = 200
      )
    }
    
    if (input$selectedCatalogType == "SBS192") {
      output$SBS192SpectrumPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(spect),
        width = 800, height = 200
      )
      output$SBS192AttributionPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(reconstructed.catalog1),
        width = 800, height = 200
      )
    }

    if (input$selectedCatalogType == "DBS78") {
      output$DBS78SpectrumPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(spect),
        width = 800, height = 200
      )
      output$DBS78AttributionPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(reconstructed.catalog1),
        width = 800, height = 200
      )
    }

    if (input$selectedCatalogType == "ID") {
      output$IDSpectrumPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(spect),
        width = 800, height = 200
      )
      output$IDAttributionPlot <- renderPlot(
        expr = ICAMS::PlotCatalog(reconstructed.catalog1),
        width = 800, height = 200
      )
    }


  })

  # When user clicks the "Remove notifications" button, all the previous
  # notifications(error, warning or message) will be removed
  observeEvent(input$remove, {
    RemoveAllNotifications(ids)
  })



}
