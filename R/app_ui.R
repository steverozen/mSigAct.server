jscode <- "
shinyjs.init = function() {
  $('#panels li a[data-value=showSpectraTab]').hide();
  $('#panels li a[data-value=sigAttributionTab]').hide();
  $('#panels li a[data-value=sigAttributionTab2]').hide();
  $('#panels li a[data-value=attributionResultsTab]').hide();
}"

#' @import shiny
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text = jscode, functions = c()),
    navbarPage(title = "mSigAct", id = "panels",
      tabPanel(title = "Home", HomeUI()),
      tabPanel(title = "Generate spectrum catalogs from VCFs",
               UploadVCFUI(), value = "generateCatalogTab"),
      tabPanel(title = "Signature attribution", UploadSpectraUI(), 
               value = "uploadSpectraTab"),
      tabPanel(title = "Show spectra", ShowSpectraUI(), 
               value = "showSpectraTab"),
      tabPanel(title = "Signature attributions", SignatureAttributionUI(), 
               value = "sigAttributionTab"),
      tabPanel(title = "Signature attributions", SignatureAttributionUI2(), 
               value = "sigAttributionTab2"),
      tabPanel(title = "Results", AttributionResultsUI(),
               value = "attributionResultsTab"),
      tabPanel(title = "Tutorial", TutorialUI(),
               value = "tutorialTab"),
      position = "fixed-top"),
    
    # Add padding because navbar pinned at the top
    tags$style(type="text/css", "body {padding-top: 70px;}"),
  )
}

#' @import shiny
TutorialUI <- function() {
  #general.guide.path <- system.file("tutorial/top.help.md", 
  #                                  package = "mSigAct.server")
  fixedPage(
    tabsetPanel(id = "helpPages",
                tabPanel(title = "General guide",
                         includeMarkdown(path = "top.help.md")),
                tabPanel(title = "Guide to generating catalogs"),
                tabPanel(title = "Guide to signature attribution"))
  )
}

#' @import shiny
AttributionResultsUI <- function() {
  fixedPage(
    tabsetPanel(id = "tabSetPanelresults",
      tabPanel(title = "Attribution counts", 
               value = "attributionCountsBest",
               DT::dataTableOutput(outputId = "exposureTable"),
               DT::dataTableOutput(outputId = "exposureTableVCF")),
      tabPanel(title = "Attribution plot", 
               value = "attributionPlotBest",
               uiOutput(outputId = "pdfview"),
               uiOutput(outputId = "pdfviewVCF"))),
  )
}

#' @import shiny
HomeUI <- function() {
  # List the first level UI elements here
  fixedPage(
    
    
    # Add a title on top the page
    titlePanel(title = p("mSigAct: Mutational Signature Activity",
                         style = "color: #337ab7"),
               windowTitle = paste0("mSigAct: Mutational Signature Activity")),


    # Add one line break
    # br(),

    # Add a short description of msigact
    h5("This web site has two main functions:"),

    h5(tags$ol(
      tags$li(actionLink(inputId = "linkToGenerateCatalogTab",
                   label = "Create and plot mutational spectrum \"catalogs\" from VCF* files")),
      br(),
      tags$li(actionLink(inputId = "linkToUploadSpectraTab",
                   label = paste0("Estimate which mutational signatures contributed to a ",
                 "mutational spectrum"))),
    )),


    h5(a(href = "https://tinyurl.com/rdzwnxd", "*VCF", target = "_blank"),
       " files contain one mutation per line, and are created ",
      "by variant callers such as ",
      a(href = "https://github.com/Illumina/strelka", "Strelka", target = "_blank"), 
      " or ",
      a(href = "https://github.com/broadgsa/gatk", "Mutect", target = "_blank")),

    # Add a link to the PCAWG7 paper about mutational signatures
    h5("For background see ",
      a(href = "https://doi.org/10.1038/s41586-020-1943-3",
        "\"The repertoire of mutational signatures in human cancer\"", target = "_blank"),
      " and ", a(href = "https://cancer.sanger.ac.uk/cosmic/signatures",
      "COSMIC Mutational Signatures", target = "_blank")),

    # Add an overview picture about the Shiny interface
    img(src = "www/rozen-mut-sig-collage.png", width = "601", height = "430"),
  )
}

#' @import shiny
UploadVCFUI <- function() {

  sidebarLayout(
    
    sidebarPanel = sidebarPanel(
      fluidRow(
        column(6, AddReferenceGenome()),
        
        column(6, AddRegion())
      ),
      
      fluidRow(
        column(6, AddVariantCaller()),
        
        # Hidden at first, pops up if Variant caller is "other"
        column(6, conditionalPanel(
          condition = "input.variantCaller == 'unknown'",
          MergeSBSsAsDBSOption()))),
      fluidRow(
        column(6, UploadVCFFiles()),
        column(
          6, MyDownloadButton(
            outputId = "download",
            label = "Create catalogs",
            style="color: #fff;
                   background-color: #337ab7;
                   border-color: #2e6da4;padding:4px;"),
          uiOutput(outputId = "showSpectraFromVCF"),
          uiOutput(outputId = "sigAttributionFromVCF"))
      )), # end sidbarPanel
      
    mainPanel = mainPanel(
    # Add the next row of control widgets
    fixedRow(
      # Add text input for user to specify the sample names
      # representing different VCF files
      column(6, AddSampleNames()),

      # Add text input for user to specify the base filename
      # of the CSV and PDF files generated
      # column(6, AddBaseFilename())
     ),

    # Add the next row of control widgets
    fixedRow(
      # Add text input for user to specify the zip file name
      column(6, AddZipfileName())
      ),


    # Add a download button for user to download VCF files to test
    fixedRow(column(6,
                    downloadButton(outputId = "downloadsampleVCFs",
                                   label = "Download example VCFs"),
                    offset = 6)),
    br(),

    # Add a button for user to run analysis on example Strelka SBS VCFs
    fixedRow(column(6,
                    MyDownloadButton(outputId = "runstrelkasbsvcfs",
                                     label =
                                       paste0("Example analysis on  ",
                                              "Strelka VCFs")),
                    offset = 6)),

    # Add one line break
    br(),

    # Add a button for user to run analysis on example Mutect VCFs
    fixedRow(column(6,
                    MyDownloadButton(outputId = "runmutectvcfs",
                                     label =
                                       paste0("Example analysis on ",
                                              "Mutect VCFs")),
                    offset = 6)),

    verbatimTextOutput(outputId = "testoutput")
  ))
}

#' @import shiny
UploadSpectraUI <- function() {
  # List the first level UI elements here
  fixedPage(
    # Add the first row of control widgets
    fixedRow(
      # Add radio buttons for user to specify the reference genome
      column(6, AddReferenceGenome2()),

      # Add radio buttons for user to specify the genomic region
      # from where the catalogs were generated
      column(6, AddRegion2())
    ),

    # Add the next row of control widgets
    fixedRow(
      # Add a file upload control for user to upload spectra file
      column(6,
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
                                      label = "ID"))
             )
    ),
    
    fixedRow(column(6, offset = 6,
                    UploadSpectra(),
                    downloadButton(outputId = "downloadSampleSpectra",
                                   label = "Download example spectra"),
    )),
    
    br(),

    # Add a download button for user to download sample spectra to test
    fixedRow(column(6, offset = 6,
                    splitLayout(cellWidths = c("30%", "70%"), 
                                uiOutput(outputId = "showSpectraFromCatalog"), 
                                uiOutput(outputId = "sigAttributionFromCatalog"))
             )),
    
    br(),
    
    fixedRow(column(6, offset = 6, 
                    uiOutput(outputId = "removeButton2")))

  )
}

#' @import shiny
ShowSpectraUI <- function() {
  fixedPage(
    shinyjs::useShinyjs(),
    sidebarLayout(

      sidebarPanel(
        uiOutput(outputId = "selectSampleFromUploadedVCF"),
        uiOutput(outputId = "selectSampleFromUploadedCatalog"),
        uiOutput(outputId = "buttonToSigAttribution")
      ),

      mainPanel(
        
        
        if (FALSE) {
          tabsetPanel(type = "tabs",
                      tabPanel("SBS96", plotOutput("SBS96plot")),
                      tabPanel("SBS192", plotOutput("SBS192plot")),
                      tabPanel("SBS1536", plotOutput("SBS1536plot")),
                      tabPanel("DBS78", plotOutput("DBS78plot")),
                      tabPanel("DBS136", plotOutput("DBS136plot")),
                      tabPanel("DBS144", plotOutput("DBS144plot")),
                      tabPanel("ID", plotOutput("IDplot"))
          )
        },
        
        uiOutput(outputId = "spectraPlotFromVCF"),
        uiOutput(outputId = "spectraPlotFromCatalog")
        
      )

    )
  )
}

#' @import shiny
SignatureAttributionUI <- function() {
  fixedPage(

    sidebarLayout(

      sidebarPanel(




          
          #actionButton(inputId = "cancel", label = "Cancel"),
          br(),
          #x7 <- uiOutput(outputId = "bookmarkButton2")
          #x9 <- uiOutput(outputId = "restoreResults")
          br(),
          x9 <- uiOutput(outputId = "downloaResults"),
      ),

      mainPanel(
        shinyjs::useShinyjs(),
        uiOutput(outputId = "attributionResults")

      )

    )
  )
}

#' @import shiny
SignatureAttributionUI2 <- function() {
  fixedPage(
    sidebarLayout(
      sidebarPanel(
        uiOutput(outputId = "selectSampleFromCatalogForAttribution2"),
        uiOutput(outputId = "selectCancerType2"),
        uiOutput(outputId = "uploadedCatalogType"),
        uiOutput(outputId = "chooseSigSubsetForSampleFromCatalog2"),
        uiOutput(outputId = "addSig"),
        br(),
        uiOutput(outputId = "chooseMoreSigs"),
        uiOutput(outputId = "analysisButton"),
        ######################################################
        uiOutput(outputId = "selectSampleFromVCFForAttribution"),
        uiOutput(outputId = "selectCancerTypeOfVCF"),
        uiOutput(outputId = "chooseCatalogType"),
        uiOutput(outputId = "chooseSigSubsetForVCF"),
        uiOutput(outputId = "addSigForVCF"),
        br(),
        uiOutput(outputId = "chooseMoreSigsForVCF"),
        uiOutput(outputId = "analysisButtonForVCF")
        ),
      
      mainPanel(
        # Must use DT::dataTableOutput instead of dataTableOutput,
        # otherwise, the data table will not show up
        DT::dataTableOutput(outputId = "mytable"),
        DT::dataTableOutput(outputId = "mytableForVCF")
      )
    )
  )
}

#' @import shiny
golem_add_external_resources <- function(){
  addResourcePath(prefix = "www", directoryPath = 
                    system.file("app/www", package = "mSigAct.server"))
  addResourcePath(prefix = "SBS96", directoryPath = 
                    system.file("app/SBS96", package = "mSigAct.server"))
  addResourcePath(prefix = "SBS192", directoryPath = 
                    system.file("app/SBS192", package = "mSigAct.server"))
  addResourcePath(prefix = "DBS78", directoryPath = 
                    system.file("app/DBS78", package = "mSigAct.server"))
  addResourcePath(prefix = "ID", directoryPath = 
                    system.file("app/ID", package = "mSigAct.server"))
  tmpdir <- tempfile()
  dir.create(tmpdir)
  addResourcePath(prefix = "results", directoryPath = tmpdir)
  
  tags$head(
    golem::activate_js(),
    # Add here all the external resources
    # If you have a custom.css in the inst/app/www
    # Or for example, you can add shinyalert::useShinyalert() here
    #tags$link(rel="stylesheet", type="text/css", href="www/custom.css")
    
    #includeCSS(path = "inst/app/www/style.css")
    
    tags$style(
      HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
      )
    )
  )
}
