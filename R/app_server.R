#' @import ICAMS
#' @import shiny
app_server <- function(input, output,session) {
  # List the first level callModules here

  files <- renderText({input$file$datapath})
  ref.genome <- reactive({input$ref.genome})
  
  
  trans.ranges <- reactive(
    if (input$ref.genome == "hg19") {
      trans.ranges.GRCh37
    }
  )
    
  if (FALSE) {
    if (input$ref.genome == "hg19") {
      trans.ranges <- trans.ranges.GRCh37
    } else if (renderPrint({input$ref.genome == "hg38"})) {
      trans.ranges <- trans.ranges.GRCh38
    } else if (renderPrint({input$ref.genome == "mm10"})) {
      trans.ranges <- trans.ranges.GRCm38
    }
  }
  
  
  if (input$vcf.type == "strelka.sbs") {
    observeEvent(input$submit, 
                 StrelkaSBSVCFFilesToZipFile(files, ref.genome, trans.ranges,
                                             input$region, input$names.of.VCFs,
                                             input$output.file, input$zipfile.name))
  }
  
  #browser()
  output$value <- renderText({input$ref.genome == "hg19"})
  
}
