#' This function extracts names of columns which contain the tumor sample
#' information in Mutect VCFs.
#'
#' @param tumor.col.names A character string containing names of columns which
#'   contain the tumor sample information in Mutect VCFs (different names
#'   separated by comma) specified by user on the browser.
#'   
#' @return A character vector containing the name of column which contains the
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
#' @return A character vector containing the sample names representing different
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
    for (i in 1:length(res$error)) {
      id.error <- showNotification(ui = "Error:", action = res$error[i], 
                                   duration = NULL, type = "error")
      id$error <- id.error
    }
  } 
  
  if (!is.null(res$warning)) {
    for (i in 1:length(res$warning)) {
      id.warning <- showNotification(ui = "Warning:", action = res$warning[i], 
                                     duration = NULL, type = "warning")
      id$warning <- c(id$warning, id.warning)
    }
  } 
  
  if (!is.null(res$message)) {
    for (i in 1:length(res$message)) {
      id.message <- showNotification(ui = "Message:", action = res$message[i], 
                                     duration = NULL, type = "message")
      id$message <- id.message
    }
  }
  return(id)
}

#' Update the existing list of notification ids.
#'
#' @param old.ids A list of notification ids to be updated.
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

#' Create a run information text file from generating zip archive from VCF
#' files.
#' 
#' @inheritParams GenerateZipFileFromMutectVCFs
#' 
#' @param nrow.data A list which contains information indicating number of data
#'   lines in the VCFs (excluding  meta-information lines and header line).
#' 
#' @importFrom  stringi stri_pad
#' 
#' @importFrom tools md5sum
#' 
#' @importFrom utils packageVersion
#' 
#' @keywords internal
AddRunInformation <- 
  function(files, vcf.names, zipfile.name, vcftype, ref.genome, 
           region, nrow.data) {
  
  run.info <- 
    file(description = file.path(tempdir(), "run-information.txt"), open = "w")
  
  # Add the header information
  header <- paste0("run-information.txt file for ", paste0(zipfile.name, ".zip"), 
                   " created on ", Sys.time())
  char.length <- nchar(header)
  writeLines(paste(rep("-", char.length), collapse = ""), run.info)
  writeLines(header, run.info)
  writeLines(paste(rep("-", char.length), collapse = ""), run.info)
  
  # Add section on purpose of ICAMS software
  writeLines("", run.info)
  writeLines("### About ICAMS ###", run.info)
  writeLines(c("Analysis and visualization of experimentally elucidated mutational",
               "signatures - the kind of analysis and visualization in Boot et al.,",
               "'In-depth characterization of the cisplatin mutational signature in",
               "human cell lines and in esophageal and liver tumors', ", 
               "Genome Research 2018, https://doi.org/10.1101/gr.230219.117.",
               "'ICAMS' stands for In-depth Characterization and Analysis of",
               "Mutational Signatures. 'ICAMS' has functions to read in variant",
               "call files (VCFs) and to collate the corresponding catalogs of",
               "mutational spectra and to analyze and plot catalogs of mutational",
               'spectra and signatures. Handles both "counts-based" and ', 
               '"density-based" catalogs of mutational spectra or signatures.'), 
             run.info)
  writeLines("", run.info)
  writeLines(c("For complete documentation of ICAMS, please refer to ",
               "https://cran.rstudio.com/web/packages/ICAMS/index.html"), run.info)
  writeLines("", run.info)
  writeLines(c("Shiny interface of ICAMS is available at ",
               "https://jnh01.shinyapps.io/icams/"), run.info)
  
  # Add ICAMS and R version used
  writeLines("", run.info)
  writeLines("### Version of the software ###", run.info)
  writeLines(paste0("ICAMS version: ", packageVersion("ICAMS")), run.info)
  writeLines(paste0("R version:     ", getRversion()), run.info)
  
  # Add input parameters specified by the user
  writeLines("", run.info)
  writeLines("### Input parameters ###", run.info)
  if (vcftype == "strelka.sbs") {
    vcftype <- "Strelka SBS VCF"
  } else if (vcftype == "strelka.id") {
    vcftype <- "Strelka ID VCF"
  } else if (vcftype == "mutect") {
    vcftype <- "Mutect VCF"
  }
  
  if (ref.genome == "hg19") {
    ref.genome <- "Human GRCh37/hg19"
  } else if (ref.genome == "hg38") {
    ref.genome <- "Human GRCh38/hg38"
  } else if (ref.genome == "mm10") {
    ref.genome <- "Mouse GRCm38/mm10"
  }
  writeLines(paste0("Type of VCF:      ", vcftype), run.info)
  writeLines(paste0("Reference genome: ", ref.genome), run.info)
  writeLines(paste0("Region:           ", region), run.info)
  
  # Add input files information
  writeLines("", run.info)
  writeLines("### Input files ###", run.info)
  max.num.of.char <- max(nchar(vcf.names))
  # Add a description of the information listed for input files
  writeLines(paste0(stringi::stri_pad("Name", width = max.num.of.char,
                                      side = "right"), "  ", 
                    "Number of data lines", "  ",
                    "MD5"), run.info)
  
  num.of.file <- length(files)
  
  nrow <- sapply(nrow.data, FUN = "[[", 1)
  for (i in 1:num.of.file) {
    writeLines(paste0(stringi::stri_pad(vcf.names[i], 
                                        width = max.num.of.char,
                                 side = "right"), "  ",
                      stringi::stri_pad(nrow[i], width = 20,
                                        side = "right"), "  ",
                      tools::md5sum(files[i])), run.info)
  }
  close(run.info)
}

#' This function generates a zip archive from Mutect VCF files.
#' 
#' @inheritParams ICAMS::MutectVCFFilesToZipFile
#' 
#' @param files Character vector of file paths to the Mutect VCF files.
#' 
#' @param vcf.names Names of VCF files uploaded by the user.
#' 
#' @param zipfile.name Name of zip archive specified by the user.
#' 
#' @import ICAMS
#' 
#' @import zip
#' 
#' @importFrom utils glob2rx
#' 
#' @keywords internal
GenerateZipFileFromMutectVCFs <- function(files,
                                          zipfile,
                                          vcf.names,
                                          zipfile.name,
                                          ref.genome, 
                                          trans.ranges = NULL, 
                                          region = "unknown", 
                                          names.of.VCFs = NULL, 
                                          tumor.col.names = NA,
                                          base.filename = "",
                                          updateProgress = NULL){
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "reading and splitting VCFs")
  }
  list <- ReadAndSplitMutectVCFs(files, names.of.VCFs, tumor.col.names)
  nrow.data <- list$nrow.data
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating SBS catalogs")
  }
  SBS.catalogs <- VCFsToSBSCatalogs(list$split.vcfs$SBS, ref.genome, 
                                    trans.ranges, region)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "generating DBS catalogs")
  }
  DBS.catalogs <- VCFsToDBSCatalogs(list$split.vcfs$DBS, ref.genome, 
                                    trans.ranges, region)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.2, detail = "generating ID catalogs")
  }
  ID.catalog <- VCFsToIDCatalogs(list$split.vcfs$ID, ref.genome, 
                                 region)[[1]]
  
  catalogs <- c(SBS.catalogs, DBS.catalogs, list(catID = ID.catalog))
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "writing catalogs to CSV files")
  }
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotting catalogs to PDF files")
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
    updateProgress(value = 0.1, detail = "generating zip archive")
  }
  
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "mutect", 
                    ref.genome, region, nrow.data)
  
  file.names <- list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf|txt"), 
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
                                  vcf.names = vcfs.info$name,
                                  zipfile.name = input$zipfile.name,
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
#' @importFrom utils glob2rx
#' 
#' @keywords internal
GenerateZipFileFromStrelkaSBSVCFs <- function(files,
                                              zipfile,
                                              vcf.names,
                                              zipfile.name,
                                              ref.genome, 
                                              trans.ranges = NULL, 
                                              region = "unknown", 
                                              names.of.VCFs = NULL, 
                                              base.filename = "",
                                              updateProgress = NULL){
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "reading and splitting VCFs")
  }
  list <- ReadAndSplitStrelkaSBSVCFs(files, names.of.VCFs)
  nrow.data <- list$nrow.data
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating SBS catalogs")
  }
  SBS.catalogs <- VCFsToSBSCatalogs(list$split.vcfs$SBS.vcfs, ref.genome, 
                                    trans.ranges, region)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "generating DBS catalogs")
  }
  DBS.catalogs <- VCFsToDBSCatalogs(list$split.vcfs$DBS.vcfs, ref.genome, 
                                    trans.ranges, region)
  
  catalogs <- c(SBS.catalogs, DBS.catalogs)
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "writing catalogs to CSV files")
  }
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".csv"))
  }
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotting catalogs to PDF files")
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
    updateProgress(value = 0.1, detail = "generating zip archive")
  }
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "strelka.sbs",
                    ref.genome, region, nrow.data)
  
  file.names <- 
    list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf|txt"), 
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
                                      vcf.names = vcfs.info$name,
                                      zipfile.name = input$zipfile.name,
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
#' @importFrom utils glob2rx
#' 
#' @keywords internal
GenerateZipFileFromStrelkaIDVCFs <- function(files,
                                             zipfile,
                                             vcf.names,
                                             zipfile.name,
                                             ref.genome, 
                                             region = "unknown", 
                                             names.of.VCFs = NULL, 
                                             base.filename = "",
                                             updateProgress = NULL){
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "reading VCFs")
  }
  list0 <- ReadStrelkaIDVCFs(files, names.of.VCFs)
  list.of.vcfs <- lapply(list0, FUN = "[[", 1)
  nrow.data <- lapply(list0, FUN = "[[", 2)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating ID catalog")
  }
  list <- VCFsToIDCatalogs(list.of.vcfs, ref.genome, region)
  
  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.4, detail = "writing catalog to CSV file")
  }
  WriteCatalog(list$catalog, file = paste0(output.file, "catID.csv"))

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotting catalog to PDF file")
  }
  PlotCatalogToPdf(list$catalog, file = paste0(output.file, "catID.pdf"))
  
  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating zip archive")
  }
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "strelka.id",
                    ref.genome, region, nrow.data)
  
  file.names <- 
    list.files(path = tempdir(), pattern = glob2rx("*.csv|pdf|txt"), 
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
  progress <- shiny::Progress$new(min = 0, max = 0.8)
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
                                     vcf.names = vcfs.info$name,
                                     zipfile.name = input$zipfile.name,
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

#' Prepare test VCFs for user to test
#' 
#' @param file Path of the file to be written.
#' 
#' @import ICAMS
#' 
#' @import zip
#' 
#' @keywords internal
PrepareTestVCFs <- function(file) {
  dir <- system.file("extdata/Strelka-SBS-vcf", package = "ICAMS")
  files <- list.files(dir, full.names = TRUE)
  zip::zipr(files, zipfile = file)
}