#' @keywords internal
CancerTypeToExposureStatData <- 
  function(tumor.cohort = "PCAWG", sig.type = "SBS96") {
    
    SampleIDToCancerType <- function(PCAWGID) {
      tumor.types <- sapply(
        PCAWGID,
        function(x) {strsplit(x, split = "::", fixed = T)[[1]][1]})
      
      return(tumor.types)
    }
    
    SplitExposureByTumorType <- function(exposure, tumor.type) {
      tt <- t(exposure)
      split.tt <- split(as.data.frame(tt), tumor.type)
      rr <- lapply(split.tt, t)
      return(rr)
    }
    
    ExposureStats1TumorType <- function(exposure) {
      count.by.sig <- rowSums(exposure) 
      present.by.sig <- count.by.sig > 0
      
      return(present.by.sig)
    }
    if (sig.type == "SBS192") {
      sig.type <- "SBS96"
    }
    exposure <- PCAWG7::exposure[[tumor.cohort]][[sig.type]]
    tumor.type <- SampleIDToCancerType(colnames(exposure))
    split.exposure <- 
      SplitExposureByTumorType(exposure, tumor.type)
    rr <- lapply(split.exposure, ExposureStats1TumorType)
    return(do.call("cbind", rr))
  }


#' Get the subset of signatures for a specified cancer type from a specified tumor cohort
#'
#' @param cancer.type Cancer type of the tumor, e.g. "Biliary-AdenoCA".
#' 
#' @param tumor.cohort The cohort of tumors to get information from package PCAWG7. Can
#' be "PCAWG", "TCGA", "other.genome", "other.exome".
#' 
#' @param sig.type Type of the signature, e.g. "SBS96".
#' 
#' @param region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @return The set of signatures found present in the specified \code{cancer.type} from
#' the specified \code{tumor.cohort}.
#' 
#' @export
CancerTypeToSigSubset <- 
  function(cancer.type, tumor.cohort = "PCAWG", sig.type = "SBS96", region = "genome") {
    
    CancerTypeToSigNames <- 
      function(cancer.type, tumor.cohort = "PCAWG", sig.type = "SBS96") {
        
        tmp <- CancerTypeToExposureStatData(tumor.cohort = tumor.cohort, 
                                            sig.type = sig.type)
        rr <- rownames(tmp[tmp[ , cancer.type], ])
        return(rr)
      }
    rname <- CancerTypeToSigNames(cancer.type, tumor.cohort = tumor.cohort,
                                  sig.type = sig.type) 
    return(PCAWG7::signature[[region]][[sig.type]][ , rname])
  }
