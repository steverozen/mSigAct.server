#' @import ICAMS
#' @import shiny
#' @import shinyFiles
app_server <- function(input, output,session) {
  # List the first level callModules here
  
  volumes <- getVolumes()
  shinyDirChoose(input, "directory", roots = volumes, session = session)
  output$dir <- renderText(parseDirPath(volumes, input$directory))
  
  observeEvent(input$submit, {
    output$download.zipfile <- 
      downloadHandler(filename = paste0(input$zipfile.name, ".zip"),
                      content = function(file) {
                        zipfile <- paste0(file, input$zipfile.name, ".zip")
                        dir <- parseDirPath(volumes, input$directory)
                        trans.ranges <- reactive({
                          if (input$ref.genome == "hg19") {
                            return(ICAMS::trans.ranges.GRCh37)
                          } else if (input$ref.genome == "hg38") {
                            return(ICAMS::trans.ranges.GRCh38)
                          } else if (input$ref.genome == "mm10") {
                            return(ICAMS::trans.ranges.GRCm38)
                          }
                        })
                        
                        names.of.VCFs <- reactive({
                          if (input$names.of.VCFs == "") {
                            return(NULL)
                          } else {
                            vector <- unlist(strsplit(input$names.of.VCFs, 
                                                      ",", fixed = TRUE))
                            return(trimws(vector))
                          }
                        })
                        
                        if (input$vcf.type == "strelka.sbs") {
                          StrelkaSBSVCFFilesToZipFile(dir,
                                                      zipfile,
                                                      input$ref.genome, 
                                                      trans.ranges(),
                                                      input$region, 
                                                      names.of.VCFs(),
                                                      input$output.file)
                        } else if (input$vcf.type == "strelka.id") {
                          StrelkaIDVCFFilesToZipFile(dir,
                                                     zipfile,
                                                     input$ref.genome,
                                                     input$region,
                                                     names.of.VCFs(),
                                                     input$output.file)
                        } else if (input$vcf.type == "mutect") {
                          tumor.col.names <- reactive({
                            if (input$tumor.col.names == "NA") {
                              return (NA)
                            } else {
                              vector1 <- unlist(strsplit(input$tumor.col.names, 
                                                        ",", fixed = TRUE))
                              return(trimws(vector1))
                            }
                          })
                          MutectVCFFilesToZipFile(dir,
                                                  zipfile,
                                                  input$ref.genome,
                                                  trans.ranges(),
                                                  input$region,
                                                  names.of.VCFs(),
                                                  tumor.col.names(),
                                                  input$output.file)
                        }
                      })
    output$download.status <- eventReactive(input$submit, "Ready to download")
  }
      )

}
