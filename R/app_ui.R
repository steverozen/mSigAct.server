#' @import shiny
app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here 
    fluidPage(
      titlePanel(paste0("In-depth Characterization and Analysis ", 
                        "of Mutational Signatures ('ICAMS')")),
      hr(),
      fluidRow(
        column(4, add_type_of_vcf()),
        column(4, add_reference_genome()),
        column(4, add_region())
      ),
      add_common_parameters(),
      fluidRow(
        column(4, add_tumor_col_names()),
        column(4, fileInput(inputId = "file", 
                            label = h5(strong("Upload VCFs of same type")), 
                            multiple = TRUE),
               fluidRow(
                 column(3, actionButton(inputId = "submit", label = "Submit")),
                 column(3, textOutput(outputId = "download.status"),
                        actionButton(inputId = "download", label = "Download zip file")),
                 )
               )
      ),
      verbatimTextOutput("value")
    )
  )
}

#' @import shiny
add_type_of_vcf <- function() {
  radioButtons(inputId = "vcf.type", label = h5(strong("Type of VCF")), 
               choices = c("Strelka SBS" = "strelka.sbs",
                           "Strelka ID" = "strelka.id",
                           "Mutect" = "mutect"),
               selected = "strelka.sbs")
}

#' @import shiny
add_reference_genome <- function(){
  radioButtons(inputId = "ref.genome", label = h5(strong("Reference genome")), 
               choices = c("Human GRCh37" = "hg19",
                           "Human GRCh38" = "hg38",
                           "Mouse GRCm38" = "mm10"),
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
            label = p(strong(h5("Names of columns containing", 
                              "tumor sample information")),
                      h5("(Only applicable to Mutect VCFs)")),
            value = "NA")
}

#' @import shiny
add_common_parameters <- function(){
  fluidRow(
    column(4, textInput(inputId = "names.of.VCFs", 
                        label = p(strong(h5("Names of VCFs")),
                                  h5("(You may leave it blank)")))),
    column(4, textInput(inputId = "output.file",
                        label = p(strong(h5("Base name of the files")),
                                  h5("(You may leave it blank)")))),
    column(4, textInput(inputId = "zipfile.name",
                        label = p(strong(h5("Name of the zip file")),
                                  h5("(You may use the default value)")),
                        value = "output"))
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
