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

#' Catch errors, warnings and messages generated when executing an R expression
#'
#' @param expr An R expression to execute.
#'
#' @return A list containing errors, warnings and messages which were generated
#'   when executing \code{expr} and the return value from \code{expr}.
#'
#' @keywords internal
CatchToList <- function(expr) {
  warning <- error <- message <- NULL
  retval <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }), warning = function(w) {
      warning <<- append(warning, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      message <<- append(message, conditionMessage(m))
    })

  error.info <- list(error = error, warning = warning, message = message)
  list(error.info = error.info, retval = retval)
}

#' @keywords internal
CheckInputsForSpectra <- function(input, catalog.path) {
  error <- NULL
  SBS192.check <- TRUE
  if (is.null(input$ref.genome2)) {
    error <- append(error, "Reference genome must be provided")
  }
  
  if (is.null(input$region2)) {
    error <- append(error, "Genomic region must be provided")
  }
  
  catalog <- ICAMS::ReadCatalog(catalog.path)
  
  if (nrow(catalog) == 192 && !(input$region2 %in% c("transcript", "exome"))) {
    error <- append(error, paste0("The genomic region for SBS192 catalog should",
                                  ' be either "transcript" or "exome"'))
    SBS192.check <- FALSE
  }
  
  return(list(error = error, SBS192.check = SBS192.check))
}

#' @keywords internal
CheckInputsForVCF <- function(input) {
  error <- NULL
  if (is.null(input$vcftype)) {
    error <- append(error, "Type of VCF files must be provided")
  }
  
  if (is.null(input$ref.genome)) {
    error <- append(error, "Reference genome must be provided")
  }
  
  
  if (is.null(input$region)) {
    error <- append(error, "Genomic region must be provided")
  }
  
  return(error)
}

#' @keywords internal
AddErrorMessage <- function(error) {
  id <- NULL
  
  if (!is.null(error)) {
    for (i in 1:length(error)) {
      id.error <- showNotification(ui = "Error:", action = error[i],
                                   duration = NULL, type = "error")
      id <- append(id, id.error)
    }
  }
  
  return(id)
}

#' @keywords internal
RemoveErrorMessage <- function(id) {
  sapply(ids, FUN = removeNotification)
}

#' Add notifications on the client browser using errors, warnings or messages
#' generated when executing an R expression
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


#' @keywords internal
CalculateNumberOfSpace <- getFromNamespace("CalculateNumberOfSpace", "ICAMS")

#' @keywords internal
AssignNumberOfAsterisks <- getFromNamespace("AssignNumberOfAsterisks", "ICAMS")

#' Create a run information text file from generating zip archive from VCF
#' files.
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
           region, mutation.loads, strand.bias.statistics) {

    run.info <-
      file(description = file.path(tempdir(), "run-information.txt"), open = "w")

    # Add the header information
    time.info <- strftime(Sys.time(), usetz = TRUE) # Get time zone information
    time.info1 <-
      gsub(pattern = "+", replacement = "UTC+", x = time.info, fixed = TRUE)
    header <- paste0("run-information.txt file for ", zipfile.name,
                     ".zip created on ", time.info1)
    char.length <- nchar(header)
    writeLines(paste(rep("-", char.length), collapse = ""), run.info)
    writeLines(header, run.info)
    writeLines(paste(rep("-", char.length), collapse = ""), run.info)

    # Add section on purpose of ICAMS software
    writeLines("", run.info)
    writeLines("--- About ICAMS ---", run.info)
    writeLines(c("Analysis and visualization of experimentally elucidated mutational",
                 "signatures - the kind of analysis and visualization in Boot et al.,",
                 '"In-depth characterization of the cisplatin mutational signature in',
                 'human cell lines and in esophageal and liver tumors", ',
                 "Genome Research 2018, https://doi.org/10.1101/gr.230219.117 and ",
                 '"Characterization of colibactin-associated mutational signature ',
                 'in an Asian oral squamous cell carcinoma and in other mucosal tumor types",',
                 'Genome Research 2020, https://doi.org/10.1101/gr.255620.119.',
                 '"ICAMS" stands for In-depth Characterization and Analysis of',
                 'Mutational Signatures. "ICAMS" has functions to read in variant',
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
                 "https://jnh01.shinyapps.io/ICAMS/"), run.info)

    # Add ICAMS and R version used
    writeLines("", run.info)
    writeLines("--- Version of the software ---", run.info)
    writeLines(paste0("ICAMS version: ", packageVersion("ICAMS")), run.info)
    writeLines(paste0("R version:     ", getRversion()), run.info)

    # Add input parameters specified by the user
    writeLines("", run.info)
    writeLines("--- Input parameters ---", run.info)
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
    writeLines("--- Input files ---", run.info)
    max.num.of.char <- max(nchar(vcf.names))
    # Add a description of the information listed for input files
    writeLines(paste0(stri_pad("Name", width = max.num.of.char,
                               side = "right"), "  ",
                      "# of data lines", "  ",
                      stri_pad("MD5", width = 32,
                               side = "right"), "  ",
                      "# of SBS", "  ",
                      "# of DBS", "  ",
                      "# of ID", "  ",
                      "# of discarded variants", "  "),
               run.info)

    num.of.file <- length(files)

    for (i in 1:num.of.file) {
      writeLines(paste0(stri_pad(vcf.names[i],
                                 width = max.num.of.char,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$total.variants[i],
                                 width = 15, side = "right"), "  ",
                        tools::md5sum(files[i]), "  ",
                        stri_pad(mutation.loads$SBS[i], width = 8,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$DBS[i], width = 8,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$ID[i], width = 7,
                                 side = "right"), "  ",
                        stri_pad(mutation.loads$discarded.variants[i],
                                 width = 22, side = "right")),
                 run.info)

    }

    if (FALSE) {
      # Add a disclaimer about discarded variants in the analysis
      writeLines("", run.info)
      writeLines(paste0("* Triplet and above base substitutions, ",
                        "complex indels, and variants with multiple alternative ",
                        "alleles are excluded from the analysis."), run.info)
    }

    # Add strand bias statistics for SBS12 plot
    if (!is.null(strand.bias.statistics)) {
      writeLines("", run.info)
      writeLines("--- Transcription strand bias statistics ---", run.info)
      list0 <- strand.bias.statistics
      num.of.sample <- length(names(list0))
      space.mat <- CalculateNumberOfSpace(list0)

      for (i in 1:num.of.sample) {
        transcribed.counts <- list0[[i]][, "transcribed"]
        untranscribed.counts <- list0[[i]][, "untranscribed"]
        q.values <- list0[[i]][, "q.values"]
        q.values.symbol <- lapply(q.values, FUN = AssignNumberOfAsterisks)
        q.values.sci <- formatC(q.values, format = "e", digits = 2)

        transcribed.info <- character(0)
        untranscribed.info <- character(0)
        header1 <- header2 <- character(0)
        mutation.class <- rownames(list0[[1]])

        for (j in 1:6) {
          header1 <- paste0(header1, stri_pad(mutation.class[j],
                                              width = space.mat[j, "space.total"],
                                              side = "both"), "|")

          header2 <-
            paste0(header2, " ",
                   stri_pad("counts",
                            width = space.mat[j, "space.counts"],
                            side = "right"), " ",
                   stri_pad("Q-value",
                            width = space.mat[j, "space.q.value"],
                            side = "right"), " ", "|")

          transcribed.info <-
            paste0(transcribed.info, " ",
                   stri_pad(transcribed.counts[j],
                            width = space.mat[j, "space.counts"],
                            side = "right"), " ",
                   stri_pad(q.values.sci[j],
                            width = space.mat[j, "space.q.value"],
                            side = "right"), " ", "|")

          untranscribed.info <-
            paste0(untranscribed.info, " ",
                   stri_pad(untranscribed.counts[j],
                            width = space.mat[j, "space.counts"],
                            side = "right"), " ",
                   stri_pad(ifelse(is.null(q.values.symbol[[j]]),
                                   "", q.values.symbol[[j]]),
                            width = space.mat[j, "space.q.value"],
                            side = "right"), " ", "|")
        }

        # Add description lines of the information listed for strand bias statistics
        writeLines(paste0(stri_pad("", width = 13), " |", header1), run.info)
        writeLines(paste0(stri_pad("Strand", width = 13, side = "right"), " |",
                          header2, "Sample name"), run.info)

        # Write the transcription strand bias statistics
        writeLines(paste0(stri_pad("transcribed", width = 13, side = "right"),
                          " |", transcribed.info, names(list0)[i]), run.info)
        writeLines(paste0(stri_pad("untranscribed", width = 13, side = "right"),
                          " |", untranscribed.info, names(list0)[i]), run.info)

        writeLines("", run.info)
      }

      # Add a description about the symbol denoting p-value
      writeLines(
        paste0("Legend: *Q<0.05, **Q<0.01, ***Q<0.001 (Benjamini-Hochberg ",
               "false discovery rates based on two-tailed binomial tests)"), run.info)

      # Add a note about direction of strand bias
      writeLines(paste0("Direction of strand bias: Fewer mutations on ",
                        "transcribed strand indicates that DNA damage occurred on ",
                        "pyrimidines,"), run.info)
      writeLines(paste0("                          Fewer mutations on ",
                        "untranscribed strand indicates that DNA damage occurred on ",
                        "purines."), run.info)
    }
    close(run.info)
  }

#' This function generates a zip archive from Mutect VCF files.
#'
#'
#' @param files Character vector of file paths to the Mutect VCF files.
#'
#' @param zipfile Pathname of the zip file to be created.
#'
#' @param vcf.names Names of VCF files uploaded by the user.
#'
#' @param zipfile.name Name of zip archive specified by the user.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link[ICAMS]{ICAMS}}.
#'
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link[BSgenome]{BSgenome}} object
#'   \enumerate{
#'     \item \code{\link[BSgenome.Hsapiens.1000genomes.hs37d5]{BSgenome.Hsapiens.1000genomes.hs37d5}}
#'     \item \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}
#'     \item \code{\link[BSgenome.Mmusculus.UCSC.mm10]{BSgenome.Mmusculus.UCSC.mm10}}
#'   }
#'   then the function will infer \code{trans.ranges} automatically. Otherwise,
#'   user will need to provide the necessary \code{trans.ranges}. Please refer to
#'   \code{\link[ICAMS]{TranscriptRanges}} for more details.
#'   If \code{is.null(trans.ranges)} do not add transcript range
#'   information.
#'
#' @param region A character string designating a genomic region;
#'  see \code{\link[ICAMS]{as.catalog}} and \code{\link[ICAMS]{ICAMS}}.
#'
#' @param names.of.VCFs Optional. Character vector of names of the VCF files.
#'   The order of names in \code{names.of.VCFs} should match the order of VCFs
#'   listed in \code{files}. If \code{NULL}(default), this function will remove
#'   all of the path up to and including the last path separator (if any) in
#'   \code{files} and file paths without extensions (and the leading dot) will be
#'   used as the names of the VCF files.
#'
#' @param tumor.col.name Optional. Character vector of column names in VCFs which contain
#'   the tumor sample information. The order of names in \code{tumor.col.names}
#'   should match the order of VCFs listed in \code{files}. If
#'   \code{tumor.col.names} is equal to \code{NA}(default), this function will
#'   use the 10th column in all the VCFs to calculate VAFs.
#'   See \code{\link[ICAMS]{GetMutectVAF}} for more details.
#'
#' @param base.filename Optional. The base name of the CSV and PDF files to be
#'   produced; multiple files will be generated, each ending in
#'   \eqn{x}\code{.csv} or \eqn{x}\code{.pdf}, where \eqn{x} indicates the type
#'   of catalog.
#'
#' @param updateProgress A callback function to update the progress indicator on
#'   the user interface.
#'
#' @import ICAMS
#'
#' @import zip
#'
#' @importFrom utils getFromNamespace
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
  split.vcfs <- ReadAndSplitMutectVCFs(files, names.of.VCFs, tumor.col.names)

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating SBS catalogs")
  }
  SBS.list <-
    ICAMS::VCFsToSBSCatalogs(list.of.SBS.vcfs = split.vcfs$SBS,
                             ref.genome = ref.genome,
                             region = region)

  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "generating DBS catalogs")
  }
  DBS.list <-
    ICAMS::VCFsToDBSCatalogs(list.of.DBS.vcfs = split.vcfs$DBS,
                             ref.genome = ref.genome,
                             region = region)

  if (is.function(updateProgress)) {
    updateProgress(value = 0.2, detail = "generating ID catalogs")
  }
  ID.list <-
    ICAMS::VCFsToIDCatalogs(list.of.vcfs = split.vcfs$ID,
                            ref.genome = ref.genome,
                            region = region)
  CombineAndReturnCatalogsForMutectVCFs <-
    getFromNamespace("CombineAndReturnCatalogsForMutectVCFs", "ICAMS")
  catalogs0 <-
    CombineAndReturnCatalogsForMutectVCFs(split.vcfs.list = split.vcfs,
                                          SBS.list = SBS.list,
                                          DBS.list = DBS.list,
                                          ID.list = ID.list)

  GetMutationLoadsFromMutectVCFs <-
    getFromNamespace("GetMutationLoadsFromMutectVCFs", "ICAMS")
  mutation.loads <- GetMutationLoadsFromMutectVCFs(catalogs0)
  strand.bias.statistics<- NULL

  # Retrieve the catalog matrix from catalogs0
  catalogs <- catalogs0
  catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL


  # Remove the ID counts catalog as it does not have abundance for
  # it to be transformed to density catalog
  catalogs$catID <- NULL

  # Transform the counts catalogs to density catalogs
  catalogs.density <- TransCountsCatalogToDensity(catalogs)

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "writing catalogs to CSV files")
  }

  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))

  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".counts.csv"))
  }

  # Write the density catalogs to CSV files
  for (name in names(catalogs.density)) {
    WriteCatalog(catalogs.density[[name]],
                 file = paste0(output.file, name, ".csv"))
  }

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotting catalogs to PDF files")
  }

  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".counts.pdf"))
    if (name == "catSBS192") {
      list <- PlotCatalogToPdf(catalogs[[name]],
                               file = paste0(output.file, "SBS12.counts.pdf"),
                               plot.SBS12 = TRUE)
      strand.bias.statistics <-
        c(strand.bias.statistics, list$strand.bias.statistics)
    }
  }

  # Plotting the density catalogs to PDFs
  for (name in names(catalogs.density)) {
    PlotCatalogToPdf(catalogs.density[[name]],
                     file = paste0(output.file, name, ".pdf"))
    if (name == "catSBS192.density") {
      list <- PlotCatalogToPdf(catalogs.density[[name]],
                               file = paste0(output.file, "SBS12.density.pdf"),
                               plot.SBS12 = TRUE)
      strand.bias.statistics <-
        c(strand.bias.statistics, list$strand.bias.statistics)
    }
  }

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating zip archive")
  }

  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "mutect",
                    ref.genome, region, mutation.loads, strand.bias.statistics)

  file.names <- list.files(path = tempdir(), pattern = "\\.(pdf|csv|txt)$",
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
  new.ids <- AddNotifications(res$error.info)

  # Update the notification ids
  return(UpdateNotificationIDs(ids, new.ids))
}

#' This function generates a zip archive from Strelka SBS VCF files.
#'
#' @inheritParams GenerateZipFileFromMutectVCFs
#'
#' @param files Character vector of file paths to the Strelka SBS VCF files.
#'
#' @import ICAMS
#'
#' @import zip
#'
#' @importFrom utils getFromNamespace
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
  split.vcfs <- ReadAndSplitStrelkaSBSVCFs(files, names.of.VCFs)

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating SBS catalogs")
  }
  SBS.list <-
    ICAMS::VCFsToSBSCatalogs(list.of.SBS.vcfs = split.vcfs$SBS.vcfs,
                            ref.genome = ref.genome,
                            region = region)

  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "generating DBS catalogs")
  }
  DBS.list <-
    ICAMS::VCFsToDBSCatalogs(list.of.DBS.vcfs = split.vcfs$DBS.vcfs,
                             ref.genome = ref.genome,
                             region = region)

  CombineAndReturnCatalogsForStrelkaSBSVCFs <-
    getFromNamespace("CombineAndReturnCatalogsForStrelkaSBSVCFs", "ICAMS")
  catalogs0 <-
    CombineAndReturnCatalogsForStrelkaSBSVCFs(split.vcfs.list = split.vcfs,
                                              SBS.list = SBS.list,
                                              DBS.list = DBS.list)
  GetMutationLoadsFromStrelkaSBSVCFs <-
    getFromNamespace("GetMutationLoadsFromStrelkaSBSVCFs", "ICAMS")
  mutation.loads <- GetMutationLoadsFromStrelkaSBSVCFs(catalogs0)
  strand.bias.statistics<- NULL

  # Retrieve the catalog matrix from catalogs0
  catalogs <- catalogs0
  catalogs$discarded.variants <- catalogs$annotated.vcfs <- NULL

  # Transform the counts catalogs to density catalogs
  catalogs.density <- TransCountsCatalogToDensity(catalogs)

  output.file <- ifelse(base.filename == "",
                        paste0(tempdir(), .Platform$file.sep),
                        file.path(tempdir(), paste0(base.filename, ".")))

  if (is.function(updateProgress)) {
    updateProgress(value = 0.3, detail = "writing catalogs to CSV files")
  }
  for (name in names(catalogs)) {
    WriteCatalog(catalogs[[name]],
                 file = paste0(output.file, name, ".counts.csv"))
  }

  # Write the density catalogs to CSV files
  for (name in names(catalogs.density)) {
    WriteCatalog(catalogs.density[[name]],
                 file = paste0(output.file, name, ".csv"))
  }

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "plotting catalogs to PDF files")
  }
  for (name in names(catalogs)) {
    PlotCatalogToPdf(catalogs[[name]],
                     file = paste0(output.file, name, ".counts.pdf"))
    if (name == "catSBS192") {
      list <- PlotCatalogToPdf(catalogs[[name]],
                               file = paste0(output.file, "SBS12.counts.pdf"),
                               plot.SBS12 = TRUE)
      strand.bias.statistics<- c(strand.bias.statistics,
                                 list$strand.bias.statistics)
    }
  }

  # Plotting the density catalogs to PDFs
  for (name in names(catalogs.density)) {
    PlotCatalogToPdf(catalogs.density[[name]],
                     file = paste0(output.file, name, ".pdf"))
    if (name == "catSBS192.density") {
      list <- PlotCatalogToPdf(catalogs.density[[name]],
                               file = paste0(output.file, "SBS12.density.pdf"),
                               plot.SBS12 = TRUE)
      strand.bias.statistics <-
        c(strand.bias.statistics, list$strand.bias.statistics)
    }
  }

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating zip archive")
  }
  AddRunInformation(files, vcf.names, zipfile.name, vcftype = "strelka.sbs",
                    ref.genome, region, mutation.loads, strand.bias.statistics)

  file.names <-
    list.files(path = tempdir(), pattern = "\\.(pdf|csv|txt)$",
               full.names = TRUE)
  zip::zipr(zipfile = zipfile, files = file.names)
  unlink(file.names)

  # Add ID catalog with zero counts
  catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID),
                  ncol = length(vcf.names))
  rownames(catID) <- ICAMS::catalog.row.order$ID
  colnames(catID) <- vcf.names
  catID <- as.catalog(object = catID, ref.genome = ref.genome,
                      region = region, catalog.type = "counts")

  catalogs$catID <- catID
  return(list(counts = catalogs, density = catalogs.density))
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
  result <- CatchToList(
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
  new.ids <- AddNotifications(result$error.info)

  # Update the notification ids
  updated.ids <- UpdateNotificationIDs(ids, new.ids)
  return(list(retval = result$retval, ids = updated.ids))
}

#' @importFrom dplyr bind_rows
#' @keywords internal
PrepareAttributionResults <- 
  function (input, output, input.catalog.type, file, plotdata) {
    output$attributionResults <- renderUI({
        tabsetPanel(
          tabPanel(title = "Attribution counts", 
                   uiOutput(outputId = "exposureTable"),
                   downloadButton(outputId = "downloadExposureTable")),
          tabPanel(title = "Attribution plot", uiOutput(outputId = "pdfview")))
    })
    
  cossim <- plotdata$cossim
  spect <- plotdata$spect
  QP.best.MAP.exp <- plotdata$QP.best.MAP.exp
  reconstructed.catalog <- plotdata$reconstructed.catalog
  sig.universe <- plotdata$sig.universe
  
  sigs.names <- QP.best.MAP.exp$sig.id
  sigs <- sig.universe[, sigs.names, drop = FALSE]
  colnames(sigs) <- 
    paste0(colnames(sigs), " (exposure = ", 
           round(QP.best.MAP.exp$QP.best.MAP.exp), ")")
  
  list.of.catalogs <- list(spect, reconstructed.catalog, sigs)
  
  output.file.path <- utils::tail(resourcePaths(), 1)
  table.file.name <- paste0("mSigAct-", colnames(spect), "-",
                            input.catalog.type, "-exposures.csv")
  pdf.file.name <- paste0("mSigAct-", colnames(spect), "-",
                          input.catalog.type, "-attribution-plot.pdf")
  pdf.file.path <- paste0(output.file.path, "/", pdf.file.name)
  table.file.path <- paste0(output.file.path, "/", table.file.name)
  
  PlotListOfCatalogsToPdf(list.of.catalogs, file = pdf.file.path)
  
  tbl1 <- data.frame(names = colnames(spect), count = colSums(spect), 
                     cosine.similarity = cossim)
  tbl2 <- data.frame(names = QP.best.MAP.exp$sig.id, 
                     count = QP.best.MAP.exp$QP.best.MAP.exp)
  tbl <- dplyr::bind_rows(tbl1, tbl2)
  
  src.file.path <- paste0("results", "/", pdf.file.name)
  output$pdfview <- renderUI({
    tags$iframe(style="height:600px; width:100%", 
                src= src.file.path)
  })
  
  if (input.catalog.type %in% c("SBS96", "SBS192")) {
    SBS.sig.names <- tbl$names[-1]
    urls <- COSMIC.v3.SBS.sig.links[SBS.sig.names, ]
    refs <- 
      paste0("<a href='",  urls, "' target='_blank'>", SBS.sig.names, "</a>")
  } else if (input.catalog.type == "DBS78") {
    DBS.sig.names <- tbl$names[-1]
    urls <- COSMIC.v3.DBS.sig.links[DBS.sig.names, ]
    refs <- 
      paste0("<a href='",  urls, "' target='_blank'>", DBS.sig.names, "</a>")
  } else if (input.catalog.type == "ID") {
    ID.sig.names <- tbl$names[-1]
    urls <- COSMIC.v3.ID.sig.links[ID.sig.names, ]
    refs <- 
      paste0("<a href='",  urls, "' target='_blank'>", ID.sig.names, "</a>")
    
  } 
  
  # Turns the names of signatures into HTML links
  tbl2 <- tbl
  tbl2$names[-1] <- refs
  
  output$exposureTable <- renderTable({
    tbl2
  }, sanitize.text.function = function(x) x, digits = 5)
  
  
  output$downloadExposureTable <- downloadHandler(
    filename = table.file.name,
    content = function(file) {
      utils::write.csv(tbl, file = file, na = "", row.names = FALSE)
    }
  )

  # Show the new attribution results
  shinyjs::show(id = "attributionResults")
  
  return(list(attribution.results = TRUE))
  
  #plotdata <- reactiveValues(spect = NULL, reconstructed.catalog = NULL,
  #                           sig.universe = NULL, QP.best.MAP.exp = NULL)
  
  #path <- system.file("extdata/mSigAct-sample-spectra.zip", 
  #                    package = "ICAMS.shiny")
  #file.copy(from = output.file, to = file)
}


#' This function generates a zip archive from Strelka ID VCF files.
#'
#' @inheritParams GenerateZipFileFromMutectVCFs
#'
#' @param files Character vector of file paths to the Strelka ID VCF files.
#'
#' @import ICAMS
#'
#' @import zip
#'
#' @importFrom utils getFromNamespace
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
  list.of.vcfs <- ReadStrelkaIDVCFs(files, names.of.VCFs)

  if (is.function(updateProgress)) {
    updateProgress(value = 0.1, detail = "generating ID catalog")
  }
  list <-
    ICAMS::VCFsToIDCatalogs(list.of.vcfs = list.of.vcfs,
                            ref.genome = ref.genome,
                            region = region)

  GetMutationLoadsFromStrelkaIDVCFs <-
    getFromNamespace("GetMutationLoadsFromStrelkaIDVCFs", "ICAMS")
  mutation.loads <- GetMutationLoadsFromStrelkaIDVCFs(list)
  strand.bias.statistics<- NULL

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
                    ref.genome, region, mutation.loads, strand.bias.statistics)

  file.names <- list.files(path = tempdir(), pattern = "\\.(pdf|csv|txt)$",
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
PrepareSampleVCFs <- function(file) {
  dir1 <- system.file("extdata/Strelka-SBS-vcf", package = "ICAMS")
  files1 <- list.files(dir1, full.names = TRUE)
  dir2 <- system.file("extdata/Strelka-ID-vcf", package = "ICAMS")
  files2 <- list.files(dir2, full.names = TRUE)
  dir3 <- system.file("extdata/Mutect-vcf", package = "ICAMS")
  files3 <- list.files(dir3, full.names = TRUE)
  files <- c(files1, files2, files3)
  zip::zipr(files, zipfile = file)
}

#' Prepare test catalogs for user to test
#'
#' @param file Path of the file to be written.
#'
#' @import zip
#'
#' @keywords internal
PrepareExampleSpectra <- function(file) {
  path <- system.file("extdata/mSigAct-example-spectra.zip", 
                     package = "ICAMS.shiny")
  file.copy(from = path, to = file)
}

#' @keywords internal
RunICAMSOnSampleStrelkaSBSVCFs <- function(output, file, ids) {
  input <- reactiveValues()
  dir <- system.file("extdata/Strelka-SBS-vcf", package = "ICAMS")
  datapath <- list.files(dir, full.names = TRUE)
  name <- tools::file_path_sans_ext(basename(datapath))
  input$vcf.files <-
    data.frame(name = name, datapath = datapath, stringsAsFactors = FALSE)
  input$names.of.VCFs <- paste(name, collapse = ", ")
  input$base.filename <- "HepG2"
  input$zipfile.name <- "mSigAct-test-run-Strelka-SBS-VCFs-output"
  input$ref.genome <- "hg19"
  input$region <- "genome"
  ProcessStrelkaSBSVCFs(input, output, file, ids)
}

#' @keywords internal
RunICAMSOnSampleMutectVCFs <- function(output, file, ids) {
  input <- reactiveValues()
  dir <- system.file("extdata/Mutect-vcf", package = "ICAMS")
  datapath <- list.files(dir, full.names = TRUE)
  name <- tools::file_path_sans_ext(basename(datapath))
  input$vcf.files <-
    data.frame(name = name, datapath = datapath, stringsAsFactors = FALSE)
  input$names.of.VCFs <- paste(name, collapse = ", ")
  input$tumor.col.names <- "NA"
  input$base.filename <- "HepG2"
  input$zipfile.name <- "mSigAct-test-run-Mutect-VCFs-output"
  input$ref.genome <- "hg19"
  input$region <- "genome"
  ProcessMutectVCFs(input, output, file, ids)
}

#' Transform a list of counts catalogs to a list of density catalogs
#'
#' @param list A list of counts catalogs.
#'
#' @return A list of density catalogs transformed from \code{list}.
#'
#' @keywords internal
TransCountsCatalogToDensity <- function(list) {
  # Create an empty list for storing the density catalogs
  list1 <- vector(mode = "list")

  for (name in names(list)) {
    name1 <- paste0(name, ".density")
    catalog <- list[[name]]
    catalog.density <-
      TransformCatalog(catalog, target.catalog.type = "density")
    list1[[name1]] <- catalog.density
  }
  return(list1)
}

# Quiets concerns of R CMD check about no visible binding for global variable
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("%...>%", ".", "ids"))
}

#' Plot List of catalogs to Pdf
#' 
#' @param list.of.catalogs List of catalogs in \code{\link{ICAMS}} format.
#'
#' @inheritParams ICAMS::PlotCatalogToPdf
#'
#' @keywords internal
PlotListOfCatalogsToPdf <- function(list.of.catalogs, 
                                    file, 
                                    plot.SBS12 = FALSE, 
                                    cex     = 0.8,
                                    grid    = TRUE, 
                                    upper   = TRUE, 
                                    xlabels = TRUE,
                                    ylim    = NULL) {
    old.par.tck.value <- graphics::par("tck")
    # Setting the width and length for A4 size plotting
    grDevices::pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
    graphics::par(tck = old.par.tck.value)
    
    num.of.catalogs <- length(list.of.catalogs)
    opar <- graphics::par(mfrow = c(8, 1), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))
    on.exit(graphics::par(opar))
    
    for (i in 1:num.of.catalogs) {
      catalog <- list.of.catalogs[[i]]
      num.of.samples <- ncol(catalog)
      
      for (j in 1:num.of.samples) {
        cat <- catalog[, j, drop = FALSE]
        PlotCatalog(cat, plot.SBS12 = plot.SBS12, cex = cex, grid = grid, 
                    upper = upper, xlabels = xlabels, ylim = ylim)
      }
      
    }
    
    grDevices::dev.off()
    invisible(list(plot.success = TRUE))
  }

