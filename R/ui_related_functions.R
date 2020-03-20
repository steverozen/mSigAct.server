#' @import shiny
#' @keywords internal
add_vcf_type <- function() {
  radioButtons(inputId = "vcftype", label = h5(strong("Type of VCF files")), 
               choiceNames = 
                 list(p(a(href = "https://github.com/Illumina/strelka", 
                          "Strelka"), " single base substitutions (SBS)"),
                      p(a(href = "https://github.com/Illumina/strelka", 
                          "Strelka"), " Indel (ID)"),
                      a(href = "https://github.com/broadgsa/gatk", "Mutect")),
               choiceValues = list("strelka.sbs", "strelka.id", "mutect"))
}

#' @import shiny
#' @keywords internal
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
#' @keywords internal
add_reference_genome <- function(){
  radioButtons(inputId = "ref.genome", 
               label = h5(strong("Reference genome")), 
               choiceNames = 
                 list("Human GRCh37/hg19", 
                      "Human GRCh38/hg38",
                      "Mouse GRCm38/mm10"),
               choiceValues = 
                 list("hg19", "hg38", "mm10"))
}

#' @import shiny
#' @keywords internal
add_region <- function(){
  radioButtons(inputId = "region", 
               label = h5(strong("Genomic region")), 
               choices = 
                 c("unknown", "genome", "exome", "transcript"))
}

#' @import shiny
#' @keywords internal
add_sample_names <- function() {
  textInput(inputId = "names.of.VCFs", 
            label = p(strong(h5("Sample names")),
                      h5("(e.g. name1, name2, name3...If blank, ", 
                         "then use the names of the VCF files.)")))
}

#' @import shiny
#' @keywords internal
add_base_filename <- function() {
  textInput(inputId = "base.filename",
            label = 
              p(strong(h5("Base name of the CSV and PDF files to create")),
                h5("(e.g. name1, or you may leave it blank)")))
}

#' @import shiny
#' @keywords internal
add_zipfile_name <- function() {
  textInput(inputId = "zipfile.name",
            label = p(strong(h5("Name of zip file to create"))),
            value = "ICAMS.output")
}

#' @import shiny
#' @keywords internal
upload_vcf_files <- function() {
  fileInput(inputId = "vcf.files", label = "Choose VCF Files", 
            multiple = TRUE)
}