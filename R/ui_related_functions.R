#' @import shiny
#' @keywords internal
AddVariantCaller <- function() {
  radioButtons(inputId = "variantCaller", 
               label = h5(strong("Variant caller"), 
                          style = "color: #337ab7"),
               choiceNames = 
                 list(a(href = "https://github.com/Illumina/strelka",
                        target = "_blank", "Strelka"),
                      a(href = "https://github.com/broadgsa/gatk", 
                        target = "_blank","Mutect"),
                      "Others"),
               choiceValues = list("strelka", "mutect", "unknown"),
               selected = character(0))
}

MergeSBSsAsDBSOption <- function() {
  radioButtons(inputId = "mergeSBS", 
               label = h5(strong("Merge adjacent SBSs as DBS"), 
                          style = "color: #337ab7"),
               choiceNames = list("yes", "no"),
               choiceValues = list("yes", "no"),
               selected = character(0))
}


#' @import shiny
#' @keywords internal
AddCatalogType <- function() {
  radioButtons(inputId = "catalogType", 
               label = h5(strong("Type of Catalogs"), 
                          style = "color: #337ab7"),
               choiceNames = 
                 list("SBS96", "SBS192", "SBS1536", "DBS78", 
                      "DBS136", "DBS144", "ID"),
               choiceValues = list("SBS96", "SBS192", "SBS1536", "DBS78", 
                                   "DBS136", "DBS144", "ID"),
               selected = character(0))
}

#' @import shiny
#' @keywords internal
AddTumorColNames <- function() {
  textInput(inputId = "tumor.col.names", 
            label = p(strong(h5("Names of columns in Mutect VCFs containing", 
                                "tumor sample information",
                                style = "color: #337ab7")),
                      h5("(e.g. colname1, colname2, colname3 ... ", 
                         "If value is NA(default), the program will use the ", 
                         "10th column in all the VCFs to calculate VAFs)")),
            value = "NA")
}

#' @import shiny
#' @keywords internal
AddReferenceGenome <- function(){
  radioButtons(inputId = "ref.genome", 
               label = h5(strong("Reference genome", style = "color: #337ab7")), 
               choiceNames = 
                 list("Human GRCh37/hg19", 
                      "Human GRCh38/hg38",
                      "Mouse GRCm38/mm10"),
               choiceValues = 
                 list("hg19", "hg38", "mm10"),
               selected = character(0))
}

#' @import shiny
#' @keywords internal
AddReferenceGenome2 <- function(){
  radioButtons(inputId = "ref.genome2", 
               label = h5(strong("Reference genome", style = "color: #337ab7")), 
               choiceNames = 
                 list("Human GRCh37/hg19", 
                      "Human GRCh38/hg38",
                      "Mouse GRCm38/mm10"),
               choiceValues = 
                 list("hg19", "hg38", "mm10"),
               selected = character(0))
}

#' @import shiny
#' @keywords internal
AddRegion <- function(){
  radioButtons(inputId = "region", 
               label = h5(strong("Genomic region"), style = "color: #337ab7"), 
               choices = 
                 # c("genome", "exome", "transcript", "unknown"),
                 # The only way to get a transcript only VCF file
                 # is if someone filters it computationally.
                 # If someone did that it would be up to them
                 # to understand the density spectrum plots for
                 # SBS96 and DBS78 are off.
                 c("genome", "exome", "unknown"),
               selected = character(0))
}

#' @import shiny
#' @keywords internal
AddRegion2 <- function(){
  radioButtons(inputId = "region2", 
               label = h5(strong("Genomic region"), style = "color: #337ab7"), 
               choices = 
                # c("genome", "exome", "transcript", "unknown"),
                # It is conceivable that someone uploads an
                # SBS96, or DBS78 spectrum that is only from a transcript,
                # e.g.because it was collapsed from a SBS192
                # spectrum, but there is no SBS96 count signature
                # derived from transcript only opportunity
                 c("genome", "exome", "unknown"),
               selected = character(0))
}

#' @import shiny
#' @keywords internal
AddSampleNames <- function() {
  textInput(inputId = "names.of.VCFs", 
            label = p(strong(h5("Sample names"), style = "color: #337ab7"),
                      h5("(e.g. name1, name2, name3...If blank, ", 
                         "then use the names of the VCF files.)")))
}

#' @import shiny
#' @keywords internal
AddBaseFilename <- function() {
  textInput(inputId = "base.filename",
            label = 
              p(strong(h5("Base name of the CSV and PDF files to create", 
                          style = "color: #337ab7")),
                h5("(e.g. name1, or you may leave it blank)")))
}

#' @import shiny
#' @keywords internal
AddZipfileName <- function() {
  textInput(inputId = "zipfile.name",
            label = p(strong(h5("Name of zip file to create", 
                                style = "color: #337ab7"))),
            value = "mSigAct-catalogs-output")
}


#' @import shiny
#' @keywords internal
UploadVCFFiles <- function() {
  div(
    fileInput(inputId = "vcf.files", label = "Choose VCF Files", 
              multiple = TRUE, accept = c(".vcf", ".VCF")), style="color: #337ab7;")
}

#' @import shiny
#' @keywords internal
UploadSpectra <- function() {
  div(
    fileInput(inputId = "upload.spectra", label = "Choose spectra file",
              accept = c(".csv", ".CSV")), style="color: #337ab7;")
}

#' @import shiny
#' @keywords internal
LoadExampleSpectra <- function() {
  textOutput(outputId = "exampleSpectraText")
  tags$head(tags$style("#exampleSpectraText{color: #337ab7;"))
  
  splitLayout(cellWidths = c("25%", "25%", "25%", "25%"),
              actionButton(inputId = "preloadSBS96Spectra", label = "SBS96"),
              actionButton(inputId = "preloadSBS192Spectra", label = "SBS192"),
              actionButton(inputId = "preloadDBS78Spectra", label = "DBS78"),
              actionButton(inputId = "preloadIDSpectra", label = "ID"))
}

#' @import shiny
#' @keywords internal
MyDownloadButton <- function (outputId, label = "Download", 
                              class = NULL, ...) {
  {
    aTag <- tags$a(id = outputId, 
                   class = paste("btn btn-default shiny-download-link", class), 
                   href = "", target = "_blank", download = NA, 
                   icon(NULL), label, ...)
  }
}