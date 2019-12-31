#' @import ICAMS
#' @import shiny
#' @import shinyFiles
app_server <- function(input, output,session) {
  # List the first level callModules here
  
  volumes <- getVolumes()
  shinyDirChoose(input, "directory", roots = volumes, session = session)
  output$dir <- renderText(parseDirPath(volumes, input$directory))
  
  output$showbutton <- eventReactive(input$submit, "yes")
  outputOptions(output, "showbutton", suspendWhenHidden = FALSE) 
  
  observeEvent(
    input$submit,
    output$download.zipfile <- 
      downloadHandler(filename = paste0(input$zipfile.name, ".zip"),
                      content = function(file) {
                        dir <- parseDirPath(volumes, input$directory)
                        trans.ranges <- reactive(GetTransRanges(input$ref.genome))
                        
                        names.of.VCFs <- reactive({
                          if (input$names.of.VCFs == "") {
                            return(NULL)
                          } else {
                            vector <- unlist(strsplit(input$names.of.VCFs, 
                                                      ",", fixed = TRUE))
                            return(trimws(vector))
                          }
                        })
                        if (input$vcftype == "strelka.sbs") {
                          withProgress(message = 'Making plot', value = 0, {
                            res <- 
                              CatchToList(StrelkaSBSVCFFilesToZipFile(dir,
                                                                      file,
                                                                      input$ref.genome, 
                                                                      trans.ranges(),
                                                                      input$region, 
                                                                      names.of.VCFs(),
                                                                      input$output.file))
                            AddMessage(output, res)
                            incProgress(amount = 1.0, detail = paste("Finishing"))
                          }
                            )
                          
                        } else if (input$vcftype == "strelka.id") {
                          
                          res <- CatchToList(
                            StrelkaIDVCFFilesToZipFile(dir,
                                                       file,
                                                       input$ref.genome,
                                                       input$region,
                                                       names.of.VCFs(),
                                                       input$output.file)
                          )
                          AddMessage(output, res)
                        } else if (input$vcftype == "mutect") {
                          tumor.col.names <- reactive({
                            if (input$tumor.col.names == "NA") {
                              return (NA)
                            } else {
                              vector1 <- unlist(strsplit(input$tumor.col.names, 
                                                         ",", fixed = TRUE))
                              return(trimws(vector1))
                            }
                          })
                          res <- CatchToList(
                            MutectVCFFilesToZipFile(dir,
                                                    file,
                                                    input$ref.genome,
                                                    trans.ranges(),
                                                    input$region,
                                                    names.of.VCFs(),
                                                    tumor.col.names(),
                                                    input$output.file)
                          )
                          AddMessage(output, res)
                        }
                      })
  )
}

#' @keywords internal
GetTransRanges <- function(ref.genome) {
  if (ref.genome == "hg19") {
    return(ICAMS::trans.ranges.GRCh37)
  } else if (ref.genome == "hg38") {
    return(ICAMS::trans.ranges.GRCh38)
  } else if (ref.genome == "mm10") {
    return(ICAMS::trans.ranges.GRCm38)
  }
}


#' @keywords internal
CatchToList <- function(expr) {
  warning <- error <- message <- NULL
  res <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }), warning = function(w) {
      warning <<- append(warning, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      message <<- append(message, conditionMessage(m))
    })
  rm(res)
  list(error = error, warning = warning, message = message)
}

#' @keywords internal
AddMessage <- function(output, res) {
  if (is.null(res$error)) {
    output$error <- NULL
  } else {
    output$error <- renderText(paste0("Error: ", res$error)) 
  }
  
  if (is.null(res$warning)) {
    output$warning <- NULL
  } else {
    output$warning <- renderText(paste0("Warning: ", res$warning))
  }
  
  if (is.null(res$message)) {
    output$message <- NULL
  } else {
    output$message <- renderText(res$message)
  }
  
}
