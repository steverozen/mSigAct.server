#' @import shiny
#' @import shinyFiles
#' @import shinythemes
app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here 
    fluidPage(
      theme = shinytheme("cerulean"),
      titlePanel(p(a("ICAMS", href = "https://github.com/steverozen/ICAMS"),
                   ": In-depth Characterization and Analysis of ", 
                   "Mutational Signatures")),
      
      hr(),
      fluidRow(
        column(3, add_type_of_vcf()),
        column(3, add_reference_genome()),
        column(3, add_region()),
        column(3, conditionalPanel(
          condition = "input.vcftype == 'mutect'",
          add_tumor_col_names()))
      ),
      add_common_parameters(),
      h5(strong("Please select a folder which contains the VCF files")),
      fluidRow(
        column(1,
               shinyDirButton(id = "directory", label = "Folder select", 
                              title = "Please select a folder which contains the VCF files"),
               textOutput("dir")),
        column(3,
               actionButton(inputId = "submit", label = "Update argument and submit"),
               textOutput("download.status"))
      ),
      rep_br(2),
      downloadButton(outputId = "download.zipfile", label = "Download Zip file"),
      rep_br(2),
      hr(),
      print(p("For complete documentation of ICAMS, please refer to ",
              a(href = "https://cran.rstudio.com/web/packages/ICAMS/index.html",
                "https://cran.rstudio.com/web/packages/ICAMS/index.html")))
     
    )
  )
}

#' @import shiny
add_type_of_vcf <- function() {
  radioButtons(inputId = "vcftype", label = h5(strong("Type of VCF")), 
               choiceNames = list(p(a(href = "https://github.com/Illumina/strelka", "Strelka"),
                                    " single base substitutions (SBS)"),
                                  p(a(href = "https://github.com/Illumina/strelka", "Strelka"),
                                    " Indel (ID)"),
                                  a(href = "https://github.com/broadgsa/gatk", "Mutect")),
               choiceValues = list("strelka.sbs", "strelka.id", "mutect"),
               selected = "strelka.sbs")
}

#' @import shiny
add_reference_genome <- function(){
  radioButtons(inputId = "ref.genome", label = h5(strong("Reference genome")), 
               choiceNames = list(p("Human GRCh37/hg19"), p("Human GRCh38/hg38"),
                                  p("Mouse GRCm38/mm10")),
               choiceValues = list("hg19", "hg38", "mm10"),
               selected = "hg19")
}

#' @import shiny
add_region <- function(){
  radioButtons(inputId = "region", label = h5(strong("Genomic region")), 
               choices = c("unknown", "genome", "exome", "transcript"), 
               selected = "unknown")
}

#' @import shiny
add_tumor_col_names <- function() {
  textInput(inputId = "tumor.col.names", 
            label = p(strong(h5("Names of columns in Mutect VCFs containing", 
                              "tumor sample information")),
                      h5("(e.g. colname1, colname2, colname3 ... ", 
                         "If value is NA(default), the program will use the ", 
                         "10th column in all the VCFs to calculate VAFs)")),
            value = "NA")
}

#' @import shiny
add_common_parameters <- function(){
  fluidRow(
    column(3, textInput(inputId = "names.of.VCFs", 
                        label = p(strong(h5("Sample name")),
                                  h5("(e.g. name1, name2, name3 ...You may leave", 
                                     "it blank to use the names of the VCF files.)")))),
    column(3, textInput(inputId = "output.file",
                        label = p(strong(h5("Base name of the CSV and PDF files to create")),
                                  h5("(e.g. name1, or you may leave it blank)")))),
    column(3, textInput(inputId = "zipfile.name",
                        label = p(strong(h5("Name of zip file to create")),
                                  h5("(e.g. name1, or you may use the default value", 
                                     " ICAMS.output)")),
                        value = "ICAMS.output"))
  )
}

#' @import shiny
golem_add_external_resources <- function(){
  
  addResourcePath(
    'www', system.file('app/www', package = 'ICAMS.shiny')
  )
 
  tags$head(
    golem::activate_js(),
    golem::favicon()
    # Add here all the external resources
    # If you have a custom.css in the inst/app/www
    # Or for example, you can add shinyalert::useShinyalert() here
    #tags$link(rel="stylesheet", type="text/css", href="www/custom.css")
  )
}
