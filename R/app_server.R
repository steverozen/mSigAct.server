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
  
  # Download sample VCFs when user clicks the button
  output$downloadsampleVCFs <- downloadHandler(
    filename = function() {
      "sample-VCFs.zip"
    }, 
    content = function(file) {
      PrepareSampleVCFs(file)
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
  
  # When user clicks the submit button, then the program will try to 
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
        output$selectSampleFromUploadedVCF <- renderUI(
          {
            
            sample.names <- colnames(counts.catalog[[1]])
            radioButtons(inputId = "sampleNameFromUploadedVCF", 
                         label = "Select the sample", 
                         choices = sample.names, 
                         selected = character(0))
          }
        )
        
        
        observeEvent(input$sampleNameFromUploadedVCF, {
          output$SBS96plot <- renderPlot({
            catSBS96 <- 
              counts.catalog$catSBS96[, input$sampleNameFromUploadedVCF, drop = FALSE]
            PlotCatalog(catSBS96)
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
        
        output$selectSampleForAttribution <- renderUI(
          {
            
            sample.names <- colnames(counts.catalog[[1]])
            selectInput(inputId = "selectedSampleForAttribution", 
                        label = "Select the sample", 
                        choices = sample.names)
          }
        )
        
        output$selectcancertype <- renderUI(
          {
            cancer.types <- 
              c(colnames(CancerTypeToExposureStatData()), "Unknown")
            selectInput(inputId = "selectedcancertype", 
                        label = "Select the cancer type", 
                        choices = cancer.types)
          }
        )
        
        output$choosecatalogtype <- renderUI(
          { 
            catalog.type <- c("SBS96", "DBS78", "ID")
            selectInput(inputId = "selectedcatalogtype", 
                        label = "Select the catalog type", 
                        choices = catalog.type)
          }
        )
        
        observeEvent(input$selectedcatalogtype, {
          output$choosesigsubsect <- renderUI(
            {
              catalog.type <- input$selectedcatalogtype
              sig.universe <- colnames(PCAWG7::signature$genome[[catalog.type]])
              
              checkboxGroupInput(inputId = "selectedsigsubset", 
                                 label = "Select the subset of signatures from COSMIC", 
                                 choices = sig.universe)
            }
          )
        }
        )
        
        
      } else if (input$vcftype == "strelka.id") {
        ids <<- ProcessStrelkaIDVCFs(input, output, file, ids)
      } else if (input$vcftype == "mutect") {
        ids <<- ProcessMutectVCFs(input, output, file, ids)
      }
    })
  
    # When user submit uploaded catalog for analysis, create radio buttons for
    # user to select the sample
    observeEvent(input$upload.catalogs, {
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
                         label = "Select the sample", 
                         choices = sample.names, 
                         selected = character(0))
          }
        )  
    })
    
    # When user uploads new catalogs, remove the previous plots
    observeEvent(input$upload.catalogs, {
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
        output$DBS786plot <- renderPlot({
          PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog, 
                              drop = FALSE])})
      } else if (input$catalogType == "DBS136") {
        output$DBS136plot <- renderPlot({
          PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog, 
                              drop = FALSE])})
      } else if (input$catalogType == "DBS144") {
        output$DBS1446plot <- renderPlot({
          PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog, 
                              drop = FALSE])})
      } else if (input$catalogType == "ID") {
        output$IDplot <- renderPlot({
          PlotCatalog(catalog[, input$selectedSampleFromUploadedCatalog, 
                              drop = FALSE])})
      }
    })
    
    
  
  # When user clicks the "Remove notifications" button, all the previous
  # notifications(error, warning or message) will be removed
  observeEvent(input$remove, {
    RemoveAllNotifications(ids)
  })
  
  
  
}
