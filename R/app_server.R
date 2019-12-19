#' @import ICAMS
#' @import shiny
app_server <- function(input, output,session) {
  # List the first level callModules here

  files <- renderText({input$file$datapath})
  vcf.type <- isolate(input$vcf.type)
  ref.genome <- isolate(input$ref.genome)
  if (ref.genome == "hg19") {
    trans.ranges <- trans.ranges.GRCh37
  } else if (ref.genome == "hg38") {
    trans.ranges <- trans.ranges.GRCh38
  } else if (ref.genome == "mm10") {
    trans.ranges <- trans.ranges.GRCm38
  }
  region <- isolate(input$region)
  
  names.of.VCFs <- isolate(input$names.of.VCFs)
  if (names.of.VCFs == "") {
    names.of.VCFs <- NULL
  }
  
  output.file <- isolate(input$output.file)
  zipfile.name <- isolate(input$zipfile.name)
  
  #browser()
 
  cats <-
    observeEvent(input$file, 
                 StrelkaSBSVCFFilesToZipFile(input$file$datapath, 
                                             ref.genome, 
                                             trans.ranges,
                                             region, names.of.VCFs,
                                             output.file,
                                             zipfile.name))
  
  
  browser()
  output$value <- files
  
  
}
