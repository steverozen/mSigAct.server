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
                        if (input$vcftype == "strelka.sbs") {
                          ProcessStrelkaSBSVCFs(input, output, file, volumes)
                        } else if (input$vcftype == "strelka.id") {
                          ProcessStrelkaIDVCFs(input, output, file, volumes)
                        } else if (input$vcftype == "mutect") {
                          ProcessMutectVCFs(input, output, file, volumes)
                        }
                      })
  )
}



#' @keywords internal
ProcessStrelkaSBSVCFs <- function(input, output, file, volumes) {
  dir <- parseDirPath(volumes, input$directory)
  trans.ranges <- reactive(GetTransRanges(input$ref.genome))
  names.of.VCFs <- reactive(GetNamesOfVCFs(input$names.of.VCFs))
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  progress$set(message = "Progress", value = 0)
  # Close the progress when this reactive exits (even if there's an error)
  on.exit(progress$close())
  
  # Create a callback function to update progress. Each time this is called, it
  # will increase the progress by that value and update the detail.
  updateProgress <- function(value = NULL, detail = NULL) {
    value1 <- value + progress$getValue()
    progress$set(value = value1, detail = detail)
  }
  
  res <- CatchToList(
    ICAMS:::.StrelkaSBSVCFFilesToZipFile(dir,
                                         file,
                                         input$ref.genome, 
                                         trans.ranges(),
                                         input$region, 
                                         names.of.VCFs(),
                                         input$output.file,
                                         updateProgress))
  AddMessage(res)
}

#' @keywords internal
ProcessStrelkaIDVCFs <- function(input, output, file, volumes) {
  dir <- parseDirPath(volumes, input$directory)
  names.of.VCFs <- reactive(GetNamesOfVCFs(input$names.of.VCFs))
  
  # Create a Progress object
  progress <- shiny::Progress$new(min = 0, max = 0.4)
  progress$set(message = "Progress", value = 0)
  # Close the progress when this reactive exits (even if there's an error)
  on.exit(progress$close())
  
  # Create a callback function to update progress. Each time this is called, it
  # will increase the progress by that value and update the detail.
  updateProgress <- function(value = NULL, detail = NULL) {
    value1 <- value + progress$getValue()
    progress$set(value = value1, detail = detail)
  }
  
  res <- CatchToList(
    ICAMS:::.StrelkaIDVCFFilesToZipFile(dir,
                                        file,
                                        input$ref.genome,
                                        input$region,
                                        names.of.VCFs(),
                                        input$output.file,
                                        updateProgress)
  )

  AddMessage(res)
}

#' @keywords internal
ProcessMutectVCFs <- function(input, output, file, volumes) {
  dir <- parseDirPath(volumes, input$directory)
  trans.ranges <- reactive(GetTransRanges(input$ref.genome))
  names.of.VCFs <- reactive(GetNamesOfVCFs(input$names.of.VCFs))
  tumor.col.names <- reactive(GetTumorColNames(input$tumor.col.names))
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  progress$set(message = "Progress", value = 0)
  # Close the progress when this reactive exits (even if there's an error)
  on.exit(progress$close())
  
  # Create a callback function to update progress. Each time this is called, it
  # will increase the progress by that value and update the detail.
  updateProgress <- function(value = NULL, detail = NULL) {
    value1 <- value + progress$getValue()
    progress$set(value = value1, detail = detail)
  }
  
  res <- CatchToList(
    ICAMS:::.MutectVCFFilesToZipFile(dir,
                                     file,
                                     input$ref.genome,
                                     trans.ranges(),
                                     input$region,
                                     names.of.VCFs(),
                                     tumor.col.names(),
                                     input$output.file,
                                     updateProgress)
  )
  AddMessage(res)
}

#' @keywords internal
GetTumorColNames <- function(tumor.col.names) {
  if (tumor.col.names == "NA") {
    return (NA)
  } else {
    vector1 <- unlist(strsplit(tumor.col.names, 
                               ",", fixed = TRUE))
    return(trimws(vector1))
  }
}


#' @keywords internal
GetNamesOfVCFs <- function(names.of.VCFs) {
  if (names.of.VCFs == "") {
    return(NULL)
  } else {
    vector <- unlist(strsplit(names.of.VCFs, 
                              ",", fixed = TRUE))
    return(trimws(vector))
  }
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
AddMessage <- function(res) {
  if (!is.null(res$error)) {
    showNotification(ui = "Error:", action = res$error, duration = NULL,
                     type = "error")
  } 
  
  if (!is.null(res$warning)) {
    showNotification(ui = "Warning:", action = res$warning, duration = NULL,
                     type = "warning")
  } 
  
  if (!is.null(res$message)) {
    showNotification(ui = "Message:", action = res$message, duration = NULL,
                     type = "message")
  } 
}
