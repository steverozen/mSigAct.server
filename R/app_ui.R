#' @import shiny
app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),

    navbarPage(title = "msigact",
      tabPanel("Overview", OverviewUI()),
      tabPanel("Generate catalogs from VCFs", UploadVCFUI()),
      tabPanel("Upload spectra", UploadSpectraUI()),
      tabPanel("Show spectra", ShowSpectraUI()),
      tabPanel("Signature attributions", SignatureAttributionUI())
    )
  )
}

#' @import shiny
OverviewUI <- function() {
  # List the first level UI elements here
  fixedPage(
    # Add a title on top the page
    titlePanel(title = p("msigact: mutational signature activity delineation ",
                         "in cancer", style = "color: #337ab7"),
               windowTitle = paste0("msigact: mutational signature activity ",
                                    "delineation in cancer")),

    # Add a horizontal line
    hr(),

    # Add a short description of msigact
    p("msigact is an online tool that allows users to upload multiple ",
      a(href = "https://tinyurl.com/rdzwnxd", "VCF"),
      " (Variant Call Format) files and returns a zip archive which ",
      "contains mutation catalogs and PDF plots. ",
      "VCF files from ",
      a(href = "https://github.com/Illumina/strelka", "Strelka"), " or ",
      a(href = "https://github.com/broadgsa/gatk", "Mutect"),
      " variant caller are supported.",
      "The uploaded VCFs must ", strong("all"), " come from the ",
      strong("same"), " variant caller, reference genome and region.",
      "User can also upload catalog to show the mutational spectrum ",
      "and do signature attribution."),

    # Add a link to the PCAWG7 paper about mutational signatures
    p("For background information of mutational signatures, please refer ",
      "to this paper: ",
      a(href = "https://doi.org/10.1038/s41586-020-1943-3",
        "The repertoire of mutational signatures in human cancer"), "."),

    # Add an overview picture about the Shiny interface
    # img(src = "www/msigact-overview.PNG", width = "800", height = "219"),

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
                                     label = "Submit",
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
                                   label = "Download the sample VCFs"),
                    offset = 6)),

    # Add one line break
    br(),

    # Add a button for user to run analysis on sample Strelka SBS VCFs
    fixedRow(column(6,
                    MyDownloadButton(outputId = "runstrelkasbsvcfs",
                                     label =
                                       paste0("Run analysis on two ",
                                              "1-sample Strelka SBS VCFs")),
                    offset = 6)),

    # Add one line break
    br(),

    # Add a button for user to run analysis on sample Mutect VCFs
    fixedRow(column(6,
                    MyDownloadButton(outputId = "runmutectvcfs",
                                     label =
                                       paste0("Run analysis on two ",
                                              "1-sample Mutect VCFs")),
                    offset = 6)),
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
      column(6, UploadSpectra(),
             actionButton(inputId = "submitSpectra", label = "Next",
                          style="color: #fff;
                                       background-color: #337ab7;
                                       border-color: #2e6da4")),
    ),
    
    br(),
    
    # Add a download button for user to download sample spectra to test
    fixedRow(column(6,
                    downloadButton(outputId = "downloadSampleSpectra",
                                   label = "Download sample spectra")))

  )
}

#' @import shiny
ShowSpectraUI <- function() {
  fixedPage(

    sidebarLayout(

      sidebarPanel(
        x <- uiOutput(outputId = "selectSampleFromUploadedVCF"),
        y <- uiOutput(outputId = "selectSampleFromUploadedCatalog")
      ),

      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("SBS96", plotOutput("SBS96plot")),
                    tabPanel("SBS192", plotOutput("SBS192plot")),
                    tabPanel("SBS1536", plotOutput("SBS1536plot")),
                    tabPanel("DBS78", plotOutput("DBS78plot")),
                    tabPanel("DBS136", plotOutput("DBS136plot")),
                    tabPanel("DBS144", plotOutput("DBS144plot")),
                    tabPanel("ID", plotOutput("IDplot"))
        )
      )

    )
  )
}

#' @import shiny
#' @import shinybusy add_busy_spinner
SignatureAttributionUI <- function() {
  fixedPage(

    sidebarLayout(

      sidebarPanel(
          x <- uiOutput(outputId = "selectSampleFromVCFForAttribution"),

          y <- uiOutput(outputId = "selectSampleFromCatalogForAttribution"),

          y1 <- uiOutput(outputId = "selectcancertype"),

          y2 <- uiOutput(outputId = "choosecatalogtype"),

          x3 <- uiOutput(outputId = "chooseSigSubsetForSampleFromVCF"),
          y4 <- uiOutput(outputId = "chooseSigSubsetForSampleFromCatalog"),

          x5 <- uiOutput(outputId = "analyzeButton1"),
          y5 <- uiOutput(outputId = "analyzeButton2")
      ),

      mainPanel(
        shinybusy::add_busy_spinner(spin = "fading-circle",
                                    color = "#2e6da4", 
                                    position = c("bottom-right")),
        
        
        uiOutput("sigContributionPlot")
        
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
