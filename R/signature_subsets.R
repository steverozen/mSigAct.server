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
    exposure <- PCAWG7::exposure[[tumor.cohort]][[sig.type]]
    tumor.type <- SampleIDToCancerType(colnames(exposure))
    split.exposure <- 
      SplitExposureByTumorType(exposure, tumor.type)
    rr <- lapply(split.exposure, ExposureStats1TumorType)
    return(do.call("cbind", rr))
  }

#' @keywords internal
CancerTypeToSigSubset <- 
  function(ca.type, tumor.cohort = "PCAWG", sig.type = "SBS96", region = "genome") {
    
    CancerTypeToSigNames <- 
      function(ca.type, tumor.cohort = "PCAWG", sig.type = "SBS96") {
        
        tmp <- CancerTypeToExposureStatData(tumor.cohort = tumor.cohort, 
                                            sig.type = sig.type)
        rr <- rownames(tmp[tmp[ , ca.type], ])
        return(rr)
      }
    rname <- CancerTypeToSigNames(ca.type, tumor.cohort = tumor.cohort,
                                  sig.type = sig.type) 
    return(PCAWG7::signature[[region]][[sig.type]][ , rname])
  }
