jscode <- "
shinyjs.init = function() {
  $('#panels li a[data-value=showSpectraTab]').hide();
  $('#panels li a[data-value=sigAttributionTab]').hide();
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
      tabPanel(title = "Upload spectra", UploadSpectraUI(), 
               value = "uploadSpectraTab"),
      tabPanel(title = "Show spectra", ShowSpectraUI(), 
               value = "showSpectraTab"),
      tabPanel(title = "Signature attributions", SignatureAttributionUI(), 
               value = "sigAttributionTab"),
      position = "fixed-top"),
    
    # Add padding because navbar pinned at the top
    tags$style(type="text/css", "body {padding-top: 70px;}")
  )
}

#' @import shinyhelper
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

    if (FALSE) {
      fixedRow(column(2, shinyhelper::helper(shiny::actionButton("go", "click me!"),
                                             icon = "question-circle",
                                             colour = "grey",
                                             type = "inline",
                                             content = "ClickHelp")))
    },


    # Add an overview picture about the Shiny interface
    img(src = "www/rozen-mut-sig-collage.png", width = "601", height = "430")

  )
}

#' @import shiny
UploadVCFUI <- function() {
  # List the first level UI elements here
  fixedPage(
    # Add the first row of control widgets
    fixedRow(
      # Add radio buttons for user to specify the type of VCF files
      column(6, AddVCFType()),

      # Add a conditional panel for user to specify the column names in
      # Mutect VCFs which contain the tumor sample information (if needed)
      column(6, conditionalPanel(
        condition = "input.vcftype == 'mutect'",
        AddTumorColNames()))
    ),

    # Add the next row of control widgets
    fixedRow(
      # Add radio buttons for user to specify the reference genome
      column(6, AddReferenceGenome()),

      # Add radio buttons for user to specify the genomic region
      # from where the VCFs were generated
      column(6, AddRegion())
    ),

    # Add the next row of control widgets
    fixedRow(
      # Add text input for user to specify the sample names
      # representing different VCF files
      column(6, AddSampleNames()),

      # Add text input for user to specify the base filename
      # of the CSV and PDF files generated
      column(6, AddBaseFilename())
    ),

    # Add the next row of control widgets
    fixedRow(
      # Add text input for user to specify the zip file name
      column(6, AddZipfileName()),

      # Add a file upload control for user to upload multiple VCF files
      column(6, UploadVCFFiles())),


    # Add the next row of control widgets
    fixedRow(column(6,
                    MyDownloadButton(outputId = "download",
                                     label = "Create catalogs",
                                     style="color: #fff;
                                       background-color: #337ab7;
                                       border-color: #2e6da4"),
                    offset = 6)),

    # Add one line break
    br(),

    fixedRow(column(6,
                    actionButton(inputId = "remove",
                                 label = "Remove notifications"),
                    offset = 6)),

    # Add one line break
    br(),

    # Add a download button for user to download VCF files to test
    fixedRow(column(6,
                    downloadButton(outputId = "downloadsampleVCFs",
                                   label = "Download example VCFs"),
                    offset = 6)),

    # Add one line break
    br(),

    # Add a button for user to run analysis on example Strelka SBS VCFs
    fixedRow(column(6,
                    MyDownloadButton(outputId = "runstrelkasbsvcfs",
                                     label =
                                       paste0("Example analysis on two ",
                                              "1-sample Strelka SBS VCFs")),
                    offset = 6)),

    # Add one line break
    br(),

    # Add a button for user to run analysis on example Mutect VCFs
    fixedRow(column(6,
                    MyDownloadButton(outputId = "runmutectvcfs",
                                     label =
                                       paste0("Example analysis on two ",
                                              "1-sample Mutect VCFs")),
                    offset = 6)),

    verbatimTextOutput(outputId = "testoutput")
  )
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
      column(6, offset = 6,
             UploadSpectra(),
             downloadButton(outputId = "downloadSampleSpectra",
                            label = "Download example spectra"),
             rep_br(2),
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
    
    br(),

    # Add a download button for user to download sample spectra to test
    fixedRow(column(6, offset = 6,
                    splitLayout(cellWidths = c("30%", "70%"), 
                                uiOutput(outputId = "showSpectra"), 
                                uiOutput(outputId = "sigAttribution"))
             )),
    
    br(),
    
    fixedRow(column(6, offset = 6, 
                    uiOutput(outputId = "removeButton2")))

  )
}

#' @import shiny
ShowSpectraUI <- function() {
  fixedPage(

    sidebarLayout(

      sidebarPanel(
        x <- uiOutput(outputId = "selectSampleFromUploadedVCF"),
        y <- uiOutput(outputId = "selectSampleFromUploadedCatalog"),
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
          x <- uiOutput(outputId = "selectSampleFromVCFForAttribution"),

          y <- uiOutput(outputId = "selectSampleFromCatalogForAttribution"),

          y1 <- uiOutput(outputId = "selectCancerType"),

          y2 <- uiOutput(outputId = "choosecatalogtype"),

          x3 <- uiOutput(outputId = "chooseSigSubsetForSampleFromVCF"),
          y5 <- uiOutput(outputId = "analyzeButtonOnTop"),
          #x7 <- uiOutput(outputId = "bookmarkButton3"),
          
          br(),
          y4 <- uiOutput(outputId = "chooseSigSubsetForSampleFromCatalog"),

          x5 <- uiOutput(outputId = "analyzeButtonForVCF"),
          splitLayout(cellWidths = c("25%", "75%"),
                      uiOutput(outputId = "analyzeButton2"),
                      uiOutput(outputId = "sigTestButton2")),
          
          br(),
          uiOutput(outputId = "cancelButton"),
          #actionButton(inputId = "cancel", label = "Cancel"),
          br(),
          x8 <- uiOutput(outputId = "bookmarkButton"),
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
golem_add_external_resources <- function(){

  addResourcePath(
    'www', system.file('app/www', package = 'ICAMS.shiny')
  )

  tags$head(
    golem::activate_js()
    # Add here all the external resources
    # If you have a custom.css in the inst/app/www
    # Or for example, you can add shinyalert::useShinyalert() here
    #tags$link(rel="stylesheet", type="text/css", href="www/custom.css")
  )
}
