#' This function extracts names of columns which contain the tumor sample
#' information in Mutect VCFs.
#'
#' @param tumor.col.names A character string containing names of columns which
#'   contain the tumor sample information in Mutect VCFs (different names
#'   separated by comma) specified by user on the browser.
#'   
#' @return A character vetor containing the name of column which contains the
#'   tumor sample information in Mutect VCFs.
#'   
#' @keywords internal
GetTumorColNames <- function(tumor.col.names) {
  if (tumor.col.names == "NA") {
    return (NA)
  } else {
    vector <- unlist(strsplit(tumor.col.names, 
                               ",", fixed = TRUE))
    return(trimws(vector))
  }
}


#' This function extracts the sample names representing different VCF files
#' specified by user.
#'
#' @param names.of.VCFs A character string containing the sample names
#'   representing different VCF files (different names separated by comma)
#'   specified by user on the browser.
#'   
#' @return A character vetor containing the sample names representing different
#'   VCF files.
#'   
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

#' Catch errors, warnings and messages generated when executing an R expression.
#'
#' @param expr An R expression to execute.
#'
#' @return A list containing errors, warnings and messages which were generated
#'   when executing \code{expr}.
#'
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

#' Add notifications on the client browser using errors, warnings or messages
#' generated when executing an R expression.
#'
#' @param res A list containing errors, warnings and messages which were
#'   generated when executing an R expression. See \code{\link{CatchToList}} for
#'   more details.
#'
#' @return A list containing the notification id for error, warning and message.
#'
#' @keywords internal
AddNotifications <- function(res) {
  # Create an empty list which can be used to store notification ids later
  id <- list("error" = character(0), "warning" = character(0), 
             "message" = character(0))
  
  if (!is.null(res$error)) {
    id.error <- showNotification(ui = "Error:", action = res$error, 
                                 duration = NULL, type = "error")
    id$error <- id.error
  } 
  
  if (!is.null(res$warning)) {
    id.warning <- showNotification(ui = "Warning:", action = res$warning, 
                                   duration = NULL, type = "warning")
    id$warning <- id.warning
  } 
  
  if (!is.null(res$message)) {
    id.message <- showNotification(ui = "Message:", action = res$message, 
                                   duration = NULL, type = "message")
    id$message <- id.message
  }
  return(id)
}

#' Update the existing list of notification ids.
#'
#' @param old.ids A list of notificaion ids to be updated.
#' 
#' @param new.ids A list of new notification ids.
#' 
#' @return A list of updated notification ids.
#'
#' @keywords internal
UpdateNotificationIDs <- function(old.ids, new.ids) {
  old.ids$error <- c(old.ids$error, new.ids$error)
  old.ids$warning <- c(old.ids$warning, new.ids$warning)
  old.ids$message <- c(old.ids$message, new.ids$message)
  return(old.ids)
}

#' Remove all notifications on the client browser.
#'
#' @param ids A list of notification ids which are to be removed.
#'
#' @keywords internal
RemoveAllNotifications <- function(ids) {
  sapply(ids$error, FUN = removeNotification)
  sapply(ids$warning, FUN = removeNotification)
  sapply(ids$message, FUN = removeNotification)
}

#' This function generates a zip archive from Mutect VCF files.
#' 
#' @inheritParams ICAMS::MutectVCFFilesToZipFile
#' 
#' @param files Character vector of file paths to the Mutect VCF files.
#' 
#' @import ICAMS
#' 
#' @import zip
#' 
#' @keywords internal
GenerateZipFileFromMutectVCFs <- function(files,
                                          zipfile,
                                          ref.genome, 
                                          trans.ranges = NULL, 
                                          region = "unknown", 
                                          names.of.VCFs = NULL, 
                                          tumor.col.names = NA,
                                          base.filename = "",
                                          updateProgress = NULL){
  
  catalogs <- MutectVCFFilesToCatalog(files,
                                      ref.genome,
                                      trans.ranges,
                                      region,
                                      names.of.VCFs,
                                      tumor.col.names,
                                      updateProgress)
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "wrote catalogs to CSV files")
  }
  
  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".pdf"))
    if (name == "catSBS192") {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, "SBS12.pdf"),
                       plot.SBS12 = TRUE)
    }
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotted catalogs to PDF files")
  }
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf"), 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
}

#' This function is a wrapper function processing Mutect VCF files to
#' generate a zip archive.
#'
#' @param input A list-like object used in shiny app to store the current values
#'   of all of the widgets in the app.
#'   
#' @param output A list-like object used in shiny app that stores instructions
#'   for building the R objects in the app.
#'   
#' @param file A file path (string) of a nonexistent temp file, using which the
#'   function writes the content to that file path. See
#'   \code{\link[shiny]{downloadHandler}} for more details.
#'   
#' @param ids A list containing the existing notification ids for error, warning
#'   and message.
#'   
#' @return A list of updated notification ids for error, warning and message
#'   after running this function.
#'
#' @keywords internal
ProcessMutectVCFs <- function(input, output, file, ids) {
  # vcfs.info is a data frame that contains one row for each uploaded file, 
  # and four columns "name", "size", "type" and "datapath". 
  # "name": The filename provided by the web browser.
  # "size": The size of the uploaded data, in bytes. 
  # "type": The MIME type reported by the browser.
  # "datapath": The path to a temp file that contains the data that was uploaded.
  vcfs.info <- input$vcf.files
  
  # Get the sample names specified by user
  vcf.names <- GetNamesOfVCFs(input$names.of.VCFs)
  
  if (is.null(vcf.names)) {
    # If user didn't specify sample names, then use VCF names
    # as the sample names
    names.of.VCFs <- 
      # Get VCF file names without extension
      tools::file_path_sans_ext(vcfs.info$name) 
  } else {
    names.of.VCFs <- vcf.names
  }
  
  # Get the names of columns containing tumor sample information in Mutect VCFs
  # specified by user
  tumor.col.names <- GetTumorColNames(input$tumor.col.names)
  
  # Get the base name of the CSV and PDF files to create specified by user
  base.filename <- input$base.filename
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  progress$set(message = "Progress", value = 0)
  # Close the progress when this reactive exits (even if there's an error)
  on.exit(progress$close())
  
  # Create a callback function to update progress. Each time this is called, it
  # will increase the progress by that value and update the detail
  updateProgress <- function(value = NULL, detail = NULL) {
    value1 <- value + progress$getValue()
    progress$set(value = value1, detail = detail)
  }
  
  # Catch the errors, warnings and messages and store them in a list when
  # generating a zip archive from Mutect VCFs
  res <- CatchToList(
    GenerateZipFileFromMutectVCFs(files = vcfs.info$datapath,
                                  zipfile = file,
                                  ref.genome = input$ref.genome, 
                                  trans.ranges = NULL, 
                                  region = input$region, 
                                  names.of.VCFs = names.of.VCFs, 
                                  tumor.col.names = tumor.col.names,
                                  base.filename = base.filename,
                                  updateProgress = updateProgress)
  )
  
  # Get the new notification ids
  new.ids <- AddNotifications(res)
  
  # Update the notification ids
  return(UpdateNotificationIDs(ids, new.ids))
}

#' This function generates a zip archive from Strelka SBS VCF files.
#' 
#' @inheritParams ICAMS::MutectVCFFilesToZipFile
#' 
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#' 
#' @import ICAMS
#' 
#' @import zip
#' 
#' @keywords internal
GenerateZipFileFromStrelkaSBSVCFs <- function(files,
                                              zipfile,
                                              ref.genome, 
                                              trans.ranges = NULL, 
                                              region = "unknown", 
                                              names.of.VCFs = NULL, 
                                              base.filename = "",
                                              updateProgress = NULL){
  
  catalogs <- StrelkaSBSVCFFilesToCatalog(files,
                                          ref.genome,
                                          trans.ranges,
                                          region,
                                          names.of.VCFs,
                                          updateProgress)
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "wrote catalogs to CSV files")
  }
  
  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".pdf"))
    if (name == "catSBS192") {
      PlotCatalogToPdf(catalogs[[name]],
                       file = paste0(output.file, "SBS12.pdf"),
                       plot.SBS12 = TRUE)
    }
  }
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotted catalogs to PDF files")
  }
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf"), 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
}

#' This function is a wrapper function processing Strelka SBS VCF files to
#' generate a zip archive.
#' 
#' @inheritParams ProcessMutectVCFs
#'   
#' @return A list of updated notification ids for error, warning and message
#'   after running this function.
#'
#' @keywords internal
ProcessStrelkaSBSVCFs <- function(input, output, file, ids) {
  # vcfs.info is a data frame that contains one row for each uploaded file, 
  # and four columns "name", "size", "type" and "datapath". 
  # "name": The filename provided by the web browser.
  # "size": The size of the uploaded data, in bytes. 
  # "type": The MIME type reported by the browser.
  # "datapath": The path to a temp file that contains the data that was uploaded.
  vcfs.info <- input$vcf.files
  
  # Get the sample names specified by user
  vcf.names <- GetNamesOfVCFs(input$names.of.VCFs)
  
  if (is.null(vcf.names)) {
    # If user didn't specify sample names, then use VCF names
    # as the sample names
    names.of.VCFs <- 
      # Get VCF file names without extension
      tools::file_path_sans_ext(vcfs.info$name) 
  } else {
    names.of.VCFs <- vcf.names
  }
  
  # Get the base name of the CSV and PDF files to create specified by user
  base.filename <- input$base.filename
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  progress$set(message = "Progress", value = 0)
  # Close the progress when this reactive exits (even if there's an error)
  on.exit(progress$close())
  
  # Create a callback function to update progress. Each time this is called, it
  # will increase the progress by that value and update the detail
  updateProgress <- function(value = NULL, detail = NULL) {
    value1 <- value + progress$getValue()
    progress$set(value = value1, detail = detail)
  }
  
  # Catch the errors, warnings and messages and store them in a list when
  # generating a zip archive from Strelka SBS VCFs
  res <- CatchToList(
    GenerateZipFileFromStrelkaSBSVCFs(files = vcfs.info$datapath,
                                      zipfile = file,
                                      ref.genome = input$ref.genome, 
                                      trans.ranges = NULL, 
                                      region = input$region, 
                                      names.of.VCFs = names.of.VCFs, 
                                      base.filename = base.filename,
                                      updateProgress = updateProgress)
  )
  
  # Get the new notification ids
  new.ids <- AddNotifications(res)
  
  # Update the notification ids
  return(UpdateNotificationIDs(ids, new.ids))
}

#' This function generates a zip archive from Strelka ID VCF files.
#' 
#' @inheritParams ICAMS::MutectVCFFilesToZipFile
#' 
#' @param files Character vector of file paths to the Strelka ID VCF files.
#' 
#' @import ICAMS
#' 
#' @import zip
#' 
#' @keywords internal
GenerateZipFileFromStrelkaIDVCFs <- function(files,
                                             zipfile,
                                             ref.genome, 
                                             region = "unknown", 
                                             names.of.VCFs = NULL, 
                                             base.filename = "",
                                             updateProgress = NULL){
  
  list <- StrelkaIDVCFFilesToCatalog(files,
                                     ref.genome,
                                     region,
                                     names.of.VCFs,
                                     updateProgress)
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  WriteCatalog(list$catalog, file = paste0(output.file, "catID.csv"))
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "wrote catalogs to CSV files")
  }
  
  PlotCatalogToPdf(list$catalog, file = paste0(output.file, "catID.pdf"))
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotted catalogs to PDF files")
  }
  
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf"), 
                           full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)
}

#' This function is a wrapper function processing Strelka ID VCF files to
#' generate a zip archive.
#' 
#' @inheritParams ProcessMutectVCFs
#'   
#' @return A list of updated notification ids for error, warning and message
#'   after running this function.
#'
#' @keywords internal
ProcessStrelkaIDVCFs <- function(input, output, file, ids) {
  # vcfs.info is a data frame that contains one row for each uploaded file, 
  # and four columns "name", "size", "type" and "datapath". 
  # "name": The filename provided by the web browser.
  # "size": The size of the uploaded data, in bytes. 
  # "type": The MIME type reported by the browser.
  # "datapath": The path to a temp file that contains the data that was uploaded.
  vcfs.info <- input$vcf.files
  
  # Get the sample names specified by user
  vcf.names <- GetNamesOfVCFs(input$names.of.VCFs)
  
  if (is.null(vcf.names)) {
    # If user didn't specify sample names, then use VCF names
    # as the sample names
    names.of.VCFs <- 
      # Get VCF file names without extension
      tools::file_path_sans_ext(vcfs.info$name) 
  } else {
    names.of.VCFs <- vcf.names
  }
  
  # Get the base name of the CSV and PDF files to create specified by user
  base.filename <- input$base.filename
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  progress$set(message = "Progress", value = 0)
  # Close the progress when this reactive exits (even if there's an error)
  on.exit(progress$close())
  
  # Create a callback function to update progress. Each time this is called, it
  # will increase the progress by that value and update the detail
  updateProgress <- function(value = NULL, detail = NULL) {
    value1 <- value + progress$getValue()
    progress$set(value = value1, detail = detail)
  }
  
  # Catch the errors, warnings and messages and store them in a list when
  # generating a zip archive from Strelka ID VCFs
  res <- CatchToList(
    GenerateZipFileFromStrelkaIDVCFs(files = vcfs.info$datapath,
                                     zipfile = file,
                                     ref.genome = input$ref.genome, 
                                     region = input$region, 
                                     names.of.VCFs = names.of.VCFs, 
                                     base.filename = base.filename,
                                     updateProgress = updateProgress)
  )
  
  # Get the new notification ids
  new.ids <- AddNotifications(res)
  
  # Update the notification ids
  return(UpdateNotificationIDs(ids, new.ids))
}