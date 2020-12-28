jscode <- "
shinyjs.init = function() {
  $('#panels li a[data-value=showSpectraTab]').hide();
  $('#panels li a[data-value=sigAttributionTab]').hide();
  $('#panels li a[data-value=attributionResultsTab]').hide();
}"

#' @import shiny
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text = jscode, functions = c()),
    navbarPage(title = tags$b("mSigAct"), id = "panels",
      tabPanel(title = tags$b("Home"), 
               HomeUI()),
      tabPanel(title = tags$b("Generate spectrum", 
                              tags$br(), "catalogs from VCFs"),
               UploadVCFUI(), 
               value = "generateCatalogTab"),
      tabPanel(title = tags$b("Upload spectra for", tags$br(),
                              "signature attribution"), 
               UploadSpectraUI(), 
               value = "uploadSpectraTab"),
      tabPanel(title = tags$b("Show spectra"), ShowSpectraUI(), 
               value = "showSpectraTab"),
      tabPanel(title = tags$b("Get signature attributions"), 
               SignatureAttributionUI(), 
               value = "sigAttributionTab"),
      tabPanel(title = tags$b("Results"), 
               AttributionResultsUI(),
               value = "attributionResultsTab"),
      tabPanel(title = tags$b("Guides (help)"), 
               TutorialUI(),
               value = "tutorialTab"),
      position = "fixed-top",
      windowTitle = "mSigAct: Mutational Signature Activity"),
    
    # Add padding because navbar pinned at the top
    tags$style(type="text/css", "body {padding-top: 75px;}"),
  )
}

#' @import shiny
TutorialUI <- function() {
  #general.guide.path <- system.file("tutorial/top.help.md", 
  #                                  package = "mSigAct.server")
  fixedPage(
    tabsetPanel(
      id = "helpPages", 
      tabPanel(title = tags$b("General guide"),
               includeMarkdown(path = "inst/guides/general.guide.md")),
      tabPanel(title = tags$b("Guide to generating catalogs"),
               includeMarkdown(
                 path = "inst/guides/Guide.to.generating.catalogs.md")),
      tabPanel(title = tags$b("Guide to signature attribution"),
               includeMarkdown(
                 path = "inst/guides/Guide.to.signature.attribution.md"))
    ) # tabsetPanel
  ) # fixedPage
}

#' @import shiny
AttributionResultsUI <- function() {
  fixedPage(
    br(),
    fixedRow(splitLayout(cellWidths = c("25%", "17.5%", "20%"),
                         downloadButton(outputId = "downloadExposureTable", 
                                        label = "Download exposure counts"),
                         downloadButton(outputId = "downloadPdf",
                                        label = "Download PDF"),
                         downloadButton(outputId = "downloadMAPTable",
                                        label = "Download diagnostics"))),
    br(),
    DT::dataTableOutput(outputId = "exposureTable"),
    
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
    h4("This web site has two main functions:"),

    h4(tags$ol(
      tags$li(actionLink(
        inputId = "linkToGenerateCatalogTab",
        label = "Generate mutational spectrum catalogs from VCFs* and plot them")),
      br(),
      tags$li(actionLink(
        inputId = "linkToUploadSpectraTab",
        label = paste0("Signature attribution: ", 
                       "Estimate which mutational signatures contributed to a ",
                       "mutational spectrum"))),
    )),
    
    
    h4(a(href = "https://tinyurl.com/rdzwnxd", "*VCFs", target = "_blank"),
       " are file that contain one mutation per line, and are created ",
      "by variant callers such as ",
      a(href = "https://github.com/Illumina/strelka", "Strelka", target = "_blank"), 
      " or ",
      a(href = "https://github.com/broadgsa/gatk", "Mutect", target = "_blank")),
    
    h4("To report issues, make suggestions, and request help, ",
       "post on ",
       a(href = "https://github.com/steverozen/mSigAct.server/issues",
         "GitHub", target = "_blank"),
       " or contact ",
       a(href = "mailto:steverozen@gmail.com", "steverozen@gmail.com",
         target = "_blank")),
    
    # Impressionistic collage
    img(src = "www/rozen-mut-sig-collage.png", 
        width = "500", height = "370"),
    
    # Add a link to the PCAWG7 paper about mutational signatures
    h4("For background see ",
       a(href = "https://doi.org/10.1038/s41586-020-1943-3",
         "\"The repertoire of mutational signatures in human cancer\"", target = "_blank"),
       " and ", a(href = "https://cancer.sanger.ac.uk/cosmic/signatures",
                  "COSMIC Mutational Signatures", target = "_blank"))
    
    
  )
}

#' @import shiny
UploadVCFUI <- function() {

  sidebarLayout(
    
    sidebarPanel = sidebarPanel(
      fluidRow(
        column(6, AddReferenceGenome()),
        
        column(6, AddRegionForVCF())
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
          6, 
          uiOutput(outputId = "clickToCreateCatalogs"),
          br(),
          uiOutput(outputId = "downloadZipFile"),
          br(),
          uiOutput(outputId = "showSpectraFromVCF"),
          br(),
          uiOutput(outputId = "sigAttributionFromVCF"))
      ), width = 5), # end sidbarPanel
      
    mainPanel = mainPanel(
    # Add the next row of control widgets
    fluidRow(
      # Add text input for user to specify the sample names
      # representing different VCF files
      column(6, AddSampleNames()),
     ),

    br(),
    # Add a download button for user to download VCF files to test
    fluidRow(column(6,
                    wellPanel(
                      h5(strong("Example data and analysis", style = "color: #337ab7")),
                      br(),
                      
                      downloadButton(outputId = "downloadsampleVCFs",
                                     label = "Download example VCFs"),
                      rep_br(2),
                      
                      actionButton(inputId = "runStrelkaVCFs", 
                                   label = paste0("Example analysis on  ",
                                                  "Strelka VCFs")),
                      
                      rep_br(2),
                      actionButton(inputId = "runMutectVCFs",
                                   label = paste0("Example analysis on ",
                                            "Mutect VCFs")),
                      rep_br(1))
                    ) # end of column
             ), width = 6) # end of mainPanel
  ) # end of sidebarLayout
}

#' @import shiny
golem_add_external_resources <- function(){
  addResourcePath(prefix = "www", directoryPath = 
                    system.file("app/www", package = "mSigAct.server"))
  
  addResourcePath(prefix = "ID", directoryPath = 
                    system.file("app/ID", package = "mSigAct.server"))
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
