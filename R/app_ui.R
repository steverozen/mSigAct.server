#' @import shiny
app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    navbarPage(title = "ICAMS.shiny",
      tabPanel("Overview", OverviewUI()),
      tabPanel("Upload VCFs", UploadVCFUI()),
      tabPanel("Upload Catalogs", UploadCatalogUI()),
      tabPanel("Show spectrums", ShowSpectrumUI()),
      tabPanel("Signature attributions", SignatureAttributionUI())
    )
  )
}

#' @import shiny
OverviewUI <- function() {
  # List the first level UI elements here 
  fixedPage(
    # Add a title on top the page
    titlePanel(title = p("ICAMS: In-depth Characterization and Analysis of ",
                         "Mutational Signatures", style = "color: #337ab7"),
               windowTitle = paste0("ICAMS: In-depth Characterization and ", 
                                    "Analysis of Mutational Signatures")),
    
    # Add a horizontal line
    hr(),
    
    # Add a short description of the shiny interface of ICAMS
    p("The ", a(href = "https://shiny.rstudio.com/", "Shiny"), 
      "interface of ",
      a(href = "https://cran.r-project.org/web/packages/ICAMS/index.html", "ICAMS"), 
      " allows users to upload multiple ", 
      a(href = "https://tinyurl.com/rdzwnxd", "VCF"), 
      " (Variant Call Format) files and returns a zip archive which ",
      "contains mutation catalogs and PDF plots. ",
      "Only VCF files from ",
      a(href = "https://github.com/Illumina/strelka", "Strelka"), " or ",
      a(href = "https://github.com/broadgsa/gatk", "Mutect"), 
      " variant caller are supported.",
      "The uploaded VCFs must ", strong("all"), " come from the ", 
      strong("same"), " variant caller, reference genome and region."),
    
    # Add a link to the COSMIC website about mutational signatures
    p("For background information of mutational signatures, please refer ", 
      "to this paper: ",
      a(href = "https://doi.org/10.1038/s41586-020-1943-3", 
        "The repertoire of mutational signatures in human cancer"), "."),
    
    # Add an overview picture about the Shiny interface
    img(src = "www/ICAMS.shiny-overview-v2.JPG", width = "800", height = "471"),
    
    # Add a horizontal line
    hr(),
    
    # Add a footer on the page showing the URL of ICAMS package
    p("For complete documentation of ICAMS, please refer to ",
      a(href = "https://cran.rstudio.com/web/packages/ICAMS/index.html",
        "https://cran.rstudio.com/web/packages/ICAMS/index.html"))
    
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
    
    # Add a button for user to run ICAMS on sample Strelka SBS VCFs
    fixedRow(column(6, 
                    MyDownloadButton(outputId = "runstrelkasbsvcfs", 
                                     label = 
                                       paste0("Run ICAMS on two ", 
                                              "1-sample Strelka SBS VCFs")),
                    offset = 6)),
    
    # Add one line break
    br(),
    
    # Add a button for user to run ICAMS on sample Mutect VCFs
    fixedRow(column(6, 
                    MyDownloadButton(outputId = "runmutectvcfs", 
                                     label = 
                                       paste0("Run ICAMS on two ", 
                                              "1-sample Mutect VCFs")),
                    offset = 6)),
    
    # Add a horizontal line
    hr(),
    
    # Add a footer on the page showing the URL of ICAMS package
    p("For complete documentation of ICAMS, please refer to ",
      a(href = "https://cran.rstudio.com/web/packages/ICAMS/index.html",
        "https://cran.rstudio.com/web/packages/ICAMS/index.html"))
    
  )
}

#' @import shiny
UploadCatalogUI <- function() {
  # List the first level UI elements here 
  fixedPage(
    # Add the first row of control widgets
    fixedRow(
      # Add radio buttons for user to specify the type of catalogs
      column(6, AddCatalogType()),
      
      # Add radio buttons for user to specify the reference genome
      column(6, AddReferenceGenome())
    ),
    
    # Add the next row of control widgets
    fixedRow(
      # Add radio buttons for user to specify the genomic region
      # from where the catalogs were generated
      column(6, AddRegion()),
      
      # Add a file upload control for user to upload multiple VCF files
      column(6, UploadCatalogs())
    ),
    
    # Add one line break
    br(),
    
    # Add a button for user to run ICAMS on sample Strelka SBS VCFs
    fixedRow(column(6, 
                    MyDownloadButton(outputId = "runstrelkasbsvcfs", 
                                     label = 
                                       paste0("Run ICAMS on two ", 
                                              "1-sample Strelka SBS VCFs")),
                    offset = 6)),
    
    # Add a horizontal line
    hr(),
    
    # Add a footer on the page showing the URL of ICAMS package
    p("For complete documentation of ICAMS, please refer to ",
      a(href = "https://cran.rstudio.com/web/packages/ICAMS/index.html",
        "https://cran.rstudio.com/web/packages/ICAMS/index.html"))
    
  )
}

#' @import shiny
ShowSpectrumUI <- function() {
  fixedPage(
    
    sidebarLayout(
      
      sidebarPanel(
        x <- uiOutput(outputId = "selectSampleFromUploadedVCF")
      ),
      
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("SBS96", plotOutput("SBS96plot")),
                    tabPanel("SBS192", plotOutput("SBS192plot")),
                    tabPanel("SBS1536", plotOutput("SBS1536plot")),
                    tabPanel("DBS78", plotOutput("DBS78plot")),
                    tabPanel("DBS136", plotOutput("DBS136plot")),
                    tabPanel("DBS144", plotOutput("DBS144plot"))
        )
      )
        
    )
  )
}

#' @import shiny
SignatureAttributionUI <- function() {
  fixedPage(
    
    sidebarLayout(
      
      sidebarPanel(
        x <- uiOutput(outputId = "selectSampleForAttribution"),
        
        x1 <- uiOutput(outputId = "selectcancertype"),
        
        x2 <- uiOutput(outputId = "choosecatalogtype"),
        
        x3 <- uiOutput(outputId = "choosesigsubsect"),
      ),
      
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("SBS96", plotOutput("SBS96attributionplot")),
                    tabPanel("DBS78", plotOutput("DBS78attributionplot")),
                    tabPanel("ID", plotOutput("DBS144attributionplot"))
        )
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
    golem::activate_js(),
    golem::favicon(ico = "ICAMS")
    # Add here all the external resources
    # If you have a custom.css in the inst/app/www
    # Or for example, you can add shinyalert::useShinyalert() here
    #tags$link(rel="stylesheet", type="text/css", href="www/custom.css")
  )
}
