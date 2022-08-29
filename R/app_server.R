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
    
    preselected.sigs <- reactiveVal(NULL)
    
    selected.more.sigs <- reactiveVal(NULL)
    
    sigs.for.attribution <- reactiveVal(NULL)
    
    show.spectra.tab.existing <- reactiveVal(FALSE)
    
    sig.attribution.tab.existing <- reactiveVal(FALSE)
    
    attribution.results.tab.existing <- reactiveVal(FALSE)
    
    analysis.for.uploaded.spectra <- reactiveVal(FALSE)
      
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
      more.sigs <- setdiff(colnames(total.signatures()), preselected.sigs())
      return(more.sigs)
    }) 
    
    total.signatures <- reactive(
      {
        if (is.null(input.catalog.type())) {
          return()
        }
        
        if (input.catalog.type() == "ID") {
          # The abundance information is not available for ID spectra, so 
          # we always use the human genome signatures for attribution analysis
          total.signatures <- COSMIC.v3.hg19.genome.ID.sigs
        } else {
          # If input.region() is "unknown", use genome signatures for attribution analysis
          if (input.region() == "unknown") {
            total.signatures <-
              COSMIC.v3.sigs[[input.ref.genome()]][["genome"]][[input.catalog.type()]]
          } else {
            total.signatures <-
              COSMIC.v3.sigs[[input.ref.genome()]][[input.region()]][[input.catalog.type()]]
          }
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
      preselected.sigs(NULL)
      selected.more.sigs(NULL)
      analysis.for.uploaded.spectra(FALSE)
      
      # When the catalogs have been generated, increment rv$catalogGeneratedFlag
      rv$catalogGeneratedFlag <- rv$catalogGeneratedFlag + 1
      
      example.strelka.vcf.analyzed <<- TRUE
      example.mutect.vcf.analyzed <<- NULL
      vcf.file.analyzed <<- NULL
      
      shinyjs::hide(id = "clickToCreateCatalogs")
      shinyjs::hide(id = "spectraPlotFromVCF")
      HideThreeOptionalTabs()
      
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
      preselected.sigs(NULL)
      selected.more.sigs(NULL)
      analysis.for.uploaded.spectra(FALSE)
      
      # When the catalogs have been generated, increment rv$catalogGeneratedFlag
      rv$catalogGeneratedFlag <- rv$catalogGeneratedFlag + 1
      
      example.mutect.vcf.analyzed <<- TRUE
      example.strelka.vcf.analyzed <<- NULL
      vcf.file.analyzed <<- NULL
      
      shinyjs::hide(id = "clickToCreateCatalogs")  
      shinyjs::hide(id = "spectraPlotFromVCF")
      HideThreeOptionalTabs()
      
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
      HideThreeOptionalTabs()
    }) 
    
    # When user clicks "Create catalogs" button, generate catalogs and show
    # "Download catalogs" button
    observeEvent(input$createCatalogs, {
      # Check the arguments for generating catalogs from VCF
      errors <- CheckInputsForVCF(input)
      ids$error <<- append(ids$error, AddErrorMessage(errors))
      
      req(input$variantCaller, input$ref.genome, input$region, input$vcf.files)
      
      # Check whether the uploaded files are really VCFs
      retval1 <- ReadAndCheckVCF(input)
      if (is.null(retval1)) {
        return()
      }
      
      if (input$variantCaller == "unknown") {
        if (is.null(input$mergeSBS)) {
          return()
        }
      }
      input.ref.genome(input$ref.genome)
      input.region(input$region)
      
      preselected.sigs(NULL)
      selected.more.sigs(NULL)
      analysis.for.uploaded.spectra(FALSE)
      
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
      
      # Check and insert two tabs "sigAttributionTab" and "showSpectraTab" on
      # navbarPage
      InsertTwoTabs(sig.attribution.tab.existing(), show.spectra.tab.existing())
      sig.attribution.tab.existing(TRUE)
      show.spectra.tab.existing(TRUE)
      
      shinyjs::show(id = "selectSampleFromUploadedVCF")
      
      shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
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
      
      # Check and insert two tabs "sigAttributionTab" and "showSpectraTab" on
      # navbarPage
      InsertTwoTabs(sig.attribution.tab.existing(), show.spectra.tab.existing())
      sig.attribution.tab.existing(TRUE)
      show.spectra.tab.existing(TRUE)
      
      shinyjs::show(id = "selectSampleFromUploadedVCF")
      
      shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab")
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
      PrepareSpectraPlotFromVCF(input = input, output = output,
                                list.of.catalogs = list.of.catalogs)
    })
    
    # When user selects the sample from uploaded catalog, show
    # the sample's mutational spectrum
    observeEvent(input$selectedSampleFromUploadedCatalog, {
      PrepareSpectraPlotFromCatalog(input = input,  output = output,
                                    input.catalog.type = input.catalog.type(),
                                    catalog = catalog)
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
            out <- ReadAndCheckCatalog(input = input, catalog.path = catalog.path, 
                                       input.region = input.region())
            if (is.null(out)) {
              return()
            } else {
              catalog <<- out
            }
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
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab")
    })
    
    observeEvent(input$clickToSigAttributionForVCF, {
      shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab")
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
      analysis.for.uploaded.spectra(TRUE)
      mutation.type <<- "SBS96"
      ShowTwoButtons()
      
      HideThreeOptionalTabs()
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
      analysis.for.uploaded.spectra(TRUE)
      mutation.type <<- "SBS192"
      ShowTwoButtons()
      
      HideThreeOptionalTabs()
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
      analysis.for.uploaded.spectra(TRUE)
      mutation.type <<- "DBS78"
      ShowTwoButtons()
      
      HideThreeOptionalTabs()
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
      analysis.for.uploaded.spectra(TRUE)
      mutation.type <<- "ID"
      ShowTwoButtons()
      
      HideThreeOptionalTabs()
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
      analysis.for.uploaded.spectra(TRUE)
      
      HideThreeOptionalTabs()
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
        out <- ReadAndCheckCatalog(input = input, catalog.path = catalog.path, 
                                   input.region = input.region())
        if (is.null(out)) {
          return()
        } else {
          catalog <<- out
        }
        
        input.catalog.type(CheckCatalogType(catalog))
        mutation.type <<- CheckCatalogType(catalog)
        
        # Hide the widgets for previous uploaded VCF
        output$spectraPlotFromVCF <- NULL
        output$selectSampleFromUploadedVCF <- NULL
        shinyjs::hide(id = "spectraPlotFromVCF")
        shinyjs::hide(id = "selectSampleFromUploadedVCF")
        
        shinyjs::hide(id = "addSig")
        shinyjs::hide(id = "chooseMoreSigs")
        shinyjs::hide(id = "chooseSigSubset")
        shinyjs::hide(id = "analysisButton")
        
        # Check and insert two tabs "sigAttributionTab" and "showSpectraTab" on
        # navbarPage
        InsertTwoTabs(sig.attribution.tab.existing(), show.spectra.tab.existing())
        sig.attribution.tab.existing(TRUE)
        show.spectra.tab.existing(TRUE)
        
        shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
        
        if (input.catalog.type() %in% c("SBS96", "SBS192", "DBS78", "ID")) {
          shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
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
        
        out <- ReadAndCheckCatalog(input = input, catalog.path = catalog.path, 
                                   input.region = input.region())
        if (is.null(out)) {
          return()
        } else {
          catalog <<- out
        }
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
        
        shinyjs::hide(id = "addSig")
        shinyjs::hide(id = "chooseMoreSigs")
        shinyjs::hide(id = "chooseSigSubset")
        shinyjs::hide(id = "analysisButton")
        
        output$uploadedCatalogType <- renderUI(
          {
            p(tags$b(paste0("Mutation type: ", mutation.type)))
          }
        )
        
        # Check and insert two tabs "sigAttributionTab" and "showSpectraTab" on
        # navbarPage
        InsertTwoTabs(sig.attribution.tab.existing(), show.spectra.tab.existing())
        sig.attribution.tab.existing(TRUE)
        show.spectra.tab.existing(TRUE)
        
        shinyjs::show(selector = '#panels li a[data-value=showSpectraTab]')
        shinyjs::show(selector = '#panels li a[data-value=sigAttributionTab]')
        shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                       selected = "sigAttributionTab")
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
            out <- ReadAndCheckCatalog(input = input, catalog.path = catalog.path, 
                                       input.region = input.region())
            if (is.null(out)) {
              return()
            } else {
              catalog <<- out
            }
            
            input.catalog.type(CheckCatalogType(catalog))
            mutation.type <<- CheckCatalogType(catalog)
            sample.names <- colnames(catalog)
            selectizeInput(inputId = "selectedSampleForAttribution",
                           label = "Select sample",
                           choices = sample.names,
                           options = list(
                             placeholder = 'Please select an option below',
                             onInitialize = I('function() { this.setValue(""); }')
                           ))
          }
        )
        
        CreateSelectCancerTypeWidget(output)
        
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
      
      if (input$selectedCancerType == "") {
        return()
      }
      
      if (!input.catalog.type() %in% c("SBS96", "SBS192", "DBS78", "ID")) {
        showNotification(ui = "Error:", 
                         action = paste0("Can only do signature attribution ", 
                                         "for SBS96, SBS192, DBS78 and ID"),
                         type = "error", duration = NULL)
        return()
      }
      
      if (input$selectedSampleForAttribution == "") {
        return()
      }
      
      selected.sig.universe <- ShowPreselectedSigs(input, output,input.catalog.type())
      # Each time when we call ShowPreselectedSigs(), must remember to update
      # the reactive value preselected.sigs
      preselected.sigs(selected.sig.universe)
    
      # Show the actionButton for user to add more signatures
      output$addSig <- renderUI(
        actionButton(inputId = "addMoreSigs", label = "Add more signatures",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      )
      shinyjs::show(id = "addSig")
    }, ignoreInit = TRUE) # end of observeEvent
    
    
    # When user clicks either of the two actionButtons "Show spectra", 
    # "Signature attribution" on UploadVCFUI page, then update
    # the shiny widgets used for attribution analysis
    CheckArgumentsForVCF <- reactive({
      list(input$showSpectraOfVCF, input$sigAttributionOfVCF)
    })
    
    observeEvent(
      CheckArgumentsForVCF(),
      {
        # Hide the widgets from previous analysis
        shinyjs::hide(id = "chooseCatalogType")
        shinyjs::hide(id = "addSig")
        shinyjs::hide(id = "chooseMoreSigs")
        shinyjs::hide(id = "chooseSigSubset")
        shinyjs::hide(id = "analysisButton")
        output$selectSampleForAttribution <- renderUI(
          { 
            sample.names <- colnames(list.of.catalogs[[1]])
            selectizeInput(inputId = "selectedSampleForAttribution",
                           label = "Select sample",
                           choices = sample.names,
                           options = list(
                             placeholder = 'Please select an option below',
                             onInitialize = I('function() { this.setValue(""); }')
                           ))
          }
        )
        
        CreateSelectCancerTypeWidget(output)
        
        # Hide the UI element previously generated by uploading spectra
        output$uploadedCatalogType <- NULL
        shinyjs::hide(id = "uploadedCatalogType")
        
        # Hide the previous signature aetiology table
        output$sigAetiologyTable <- NULL
        shinyjs::hide(id = "sigAetiologyTable")
      }) # end of observeEvent
    
    # Only when user selects the sample, then show the widget to choose catalog
    # type for signature attribution
    observeEvent(input$selectedSampleForAttribution, {
      req(input$selectedSampleForAttribution)
      if (analysis.for.uploaded.spectra() == TRUE) {
        shinyjs::hide(id = "addSig")
        shinyjs::hide(id = "chooseMoreSigs")
        shinyjs::hide(id = "chooseSigSubset")
        shinyjs::hide(id = "analysisButton")
        shinyjs::hide(id = "chooseCatalogType")
        output$sigAetiologyTable <- NULL
        shinyjs::hide(id = "sigAetiologyTable")
        preselected.sigs(NULL)
        selected.more.sigs(NULL)
        
        CreateSelectCancerTypeWidget(output)
        return()
      } else {
        # Remove the values from previous sample
        input.catalog.type(NULL)
        
        shinyjs::hide(id = "addSig")
        shinyjs::hide(id = "chooseMoreSigs")
        shinyjs::hide(id = "chooseSigSubset")
        shinyjs::hide(id = "analysisButton")
        
        CreateSelectCancerTypeWidget(output)
        
        # Determine the catalog types available for attribution for selected sample
        catalog.types.for.attribution <-
          DetermineCatalogTypesForAttribution(list.of.catalogs = list.of.catalogs,
                                              sample.name = input$selectedSampleForAttribution)
        
        output$chooseCatalogType <- renderUI(
          {
            selectizeInput(inputId = "selectCatalogType", 
                           label = "Select catalog type",
                           choices = catalog.types.for.attribution,
                           options = list(
                             placeholder = 'Please select an option below',
                             onInitialize = I('function() { this.setValue(""); }')
                           ))
          }
        )
        
        shinyjs::show(id = "chooseCatalogType")
      }
    })
    
    observeEvent(input$selectCatalogType, {
      req(input$selectCatalogType)
      
      # Don't update the input.catalog.type using input$selectCatalogType
      # We are doing analysis using uploaded spectra, however input$selectCatalogType
      # is related to uploaded VCF
      if (is.null(catalog)) {
        catalog.name <- paste0("cat", input.catalog.type())
        catalog <<- list.of.catalogs[[catalog.name]]
      }
      
      selected.sig.universe <- ShowPreselectedSigs(input, output, input.catalog.type())
      # Each time when we call ShowPreselectedSigs(), must remember to update
      # the reactive value preselected.sigs
      preselected.sigs(selected.sig.universe)
      
      # Show the actionButton for user to add more signatures
      output$addSig <- renderUI(
        actionButton(inputId = "addMoreSigs", label = "Add more signatures",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      )
      
      shinyjs::show(id = "addSig")
      
    }, ignoreInit = TRUE) # end of observeEvent
    
    observeEvent(input$selectCatalogType, {
      req(input$selectCatalogType)
      
      input.catalog.type(input$selectCatalogType) 
      if (is.null(catalog)) {
        catalog.name <- paste0("cat", input.catalog.type())
        catalog <<- list.of.catalogs[[catalog.name]]
      }
      selected.sig.universe <- ShowPreselectedSigs(input, output, input.catalog.type())
      # Each time when we call ShowPreselectedSigs(), must remember to update
      # the reactive value preselected.sigs
      preselected.sigs(selected.sig.universe)
      
      # Show the actionButton for user to add more signatures
      output$addSig <- renderUI(
        actionButton(inputId = "addMoreSigs", label = "Add more signatures",
                     style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
      )
      shinyjs::show(id = "addSig")
    }, ignoreInit = TRUE) # end of observeEvent
    
    # Update the catalog used for attribution if user selects another catalog type
    # for analysis
    observeEvent(input$selectCatalogType, {
      if (input$selectCatalogType == "") {
        return()
      }
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
      shinyjs::show(id = "chooseMoreSigs")
    })
    
    # Update the reactive value when user clicks the widget
    # Must use observer() instead of observeEvent because
    # observeEvent does not react to NULL by default
    observe({
      preselected.sigs(input$preselectedSigs)
    })
    
    observe({
      selected.more.sigs(input$selectedMoreSigs)
      # Give more priority for this observer, otherwise when starting with
      # new R session, the sig etiology table will not get updated after user
      # selects more signatures
    }, priority = 1)
    
    # Create a reactive expression to determine the signatures used for
    # attribution and then show the signatures aetiology information table
    sigsForAttributionChanged <- reactive({
      list(input$selectedSampleForAttribution, input$selectCatalogType, 
           input$preselectedSigs, input$selectedMoreSigs, preselected.sigs())
    })
    
    observeEvent(input$selectCatalogType, {
      sigs.for.attribution(NULL)
      output$sigAetiologyTable <- NULL
      shinyjs::hide(id = "sigAetiologyTable")
    })
    
    observeEvent(sigsForAttributionChanged(), {
      shinyjs::hide(id = "sigAetiologyTable")
      sigs.for.attribution(c(preselected.sigs(), selected.more.sigs()))
      
      if (is.null(sigs.for.attribution())) {
        return()
      }
      
      if (is.null(input.catalog.type())) {
        return()
      }
      
      if (input$selectedSampleForAttribution == "") {
        return()
      }
      
      
      if (input.catalog.type() == "SBS96") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.SBS96.sig.links), sigs.for.attribution())
      } else if (input.catalog.type() == "SBS192") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.SBS192.sig.links), sigs.for.attribution())
      } else if (input.catalog.type() == "DBS78") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.DBS78.sig.links), sigs.for.attribution())
      } else if (input.catalog.type() == "ID") {
        sigs.in.correct.order <- 
          intersect(rownames(COSMIC.v3.ID.sig.links), sigs.for.attribution())
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
                                       search = "Search in signatures and etiologies:"
                                     )
                      )
        )
      })
      
      if (!is.null(sigs.for.attribution())) {
        shinyjs::show(id = "sigAetiologyTable")
      }
      
      # Don't show "Analyze" button if no signatures were selected in the beginning
      if (!is.null(sigs.for.attribution())) {
        output$analysisButton <- renderUI({
          actionButton(inputId = "startAnalysis", label = "Analyze",
                       style= "color: #fff; background-color: #337ab7;
                              border-color: #2e6da4; ")
        })
        shinyjs::show(id = "analysisButton")
      } # Need to use ignoreNULL here otherwise sigAetiologyTable will not change
        # if user deselect all signatures
    }, ignoreNULL = FALSE,  ignoreInit = TRUE)
    
    # Asynchronous programming starts from here
    observeEvent(input$startAnalysis, {
      if(length(sigs.for.attribution()) == 0) {
        showNotification(ui = "Error:", 
                         action = "No signatures selected for attribution analysis",
                         type = "error", duration = NULL)
        return()
      }
      
      # Don't do anything if in the middle of a run
      if(running()) {
        return(NULL)
      }
      
      # Create a Progress object
      progress <- ipc::AsyncProgress$new(session, min = 0, max = 1,
                                         message = "Analysis in progress",
                                         detail = "This may take a while...")
      
      running(TRUE)
      spect <- catalog[, input$selectedSampleForAttribution, drop = FALSE]
      catalog.type <- input.catalog.type()
      cancer.type <- input$selectedCancerType
      
      if (catalog.type == "ID") {
        sig.universe <<- 
          COSMIC.v3.hg19.genome.ID.sigs[, sigs.for.attribution(), drop = FALSE]
      } else {
        
        sig.universe <<- total.signatures()[, sigs.for.attribution(), drop = FALSE]
      }
      
      # Do the first round of cut-off if there are many signatures in the beginning
      if (ncol(sig.universe) > 15) {
        progress$inc(amount = 0.05, 
                     detail = "Trying to remove some signatures using bootstrapping")
        retval <- 
          mSigAct:::OptimizeExposureQPBootstrap(spectrum = spect,
                                                signatures = sig.universe, 
                                                mc.cores = AdjustNumberOfCores(50))
        sig.universe <<- sig.universe[, names(retval$exposure), drop = FALSE]
        progress$inc(amount = 0.05, 
                     detail = "Removed some signatures using bootstrapping")
      }
      
      sigs.prop <- mSigAct::ExposureProportions(
        mutation.type = catalog.type,
        cancer.type   = cancer.type,
        all.sigs      = sig.universe,
        must.include  = colnames(sig.universe))
      
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
          
          retval <- mSigAct:::MAPAssignActivity1(
            spect = spect,
            sigs = sig.universe,
            sigs.presence.prop = sigs.prop,
            max.level = length(sigs.prop) - 1,
            p.thresh = 0.01,
            m.opts = mSigAct::DefaultManyOpts(),
            max.mc.cores = AdjustNumberOfCores(50),
            progress.monitor = updateProgress,
            use.sparse.assign = TRUE,
            drop.low.mut.samples = FALSE
          )
          
          return(retval)
        }, seed = TRUE) %...>% {
          retval <- .
          plotdata$retval <<- retval
          
          if (is.null(retval$proposed.assignment)) {
            showNotification(ui = "Message:", 
                             action = retval$error.messages,
                             type = "message", duration = NULL)
            return()
          } else {
            MAP.best.exp <- retval$proposed.assignment
            
            reconstructed.catalog0 <- retval$proposed.reconstruction
            
            cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
            
            plotdata$cossim <<- cossim
            plotdata$spect <<- spect
            plotdata$reconstructed.catalog <<- reconstructed.catalog0
            plotdata$sig.universe <<- sig.universe
            plotdata$best.MAP.exp <<- MAP.best.exp
            retval <- 
              PrepareAttributionResults(input, output, session, 
                                        input.catalog.type(), 
                                        plotdata,
                                        attribution.results.tab.existing())
            attribution.results.tab.existing(TRUE)
            
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
    
    ########################################################################
    # Start of functions related to AttributionResultsUI
    ########################################################################
    # When user clicks the actionButton "Analyze another sample", 
    # direct user to the SignatureAttributionUI page
    observeEvent(input$analyzeMoreSample, {
      shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                     selected = "sigAttributionTab")
    })
  }, 
  error = function(err.info) {
    if (!is.null(err.info$message)) {
      showNotification(err.info$message, duration = NULL, 
                       type = "error")
    }
  }) # end of tryCatch
}
