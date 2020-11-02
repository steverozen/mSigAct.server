#' Get signature exposure for one sample with confidence interval
#'
#' @param catalog A \strong{counts} catalog as defined in \code{\link{ICAMS}}.
#'   It can only has \strong{one} column.
#'   
#' @param sig.universe The universe of signatures used to do signature attribution.
#' 
#' @param num.of.replicates The number of bootstrap replicates. 
#' 
#' @param method Method used to get the optimal solution for signature attribution.
#' 
#' @param conf.int A number specifying the required confidence interval.
#'
#' @section Value 
#' A matrix showing the signature exposure results for
#'   \code{catalog} with lower and upper bound of the sepcified confidence
#'   interval.
#'   
#' @keywords internal
GetExposureWithConfidence <- function(catalog, 
                                      sig.universe, 
                                      num.of.replicates = 1000, 
                                      method = decomposeQP,
                                      conf.int = 0.95) {
  # Get the original total mutation counts
  total.counts <- colSums(catalog)
  
  retval <- bootstrapSigExposures(m = catalog, 
                                  P = sig.universe, 
                                  R = num.of.replicates,
                                  decomposition.method = method)
  
  exposure.mean <- apply(retval$exposures, MARGIN = 1, FUN = mean)
  
  # Use function quantile to find the confidence interval of bootstrap sample values
  prob1 <- (1 - conf.int) / 2
  prob2 <- 1 - prob1
  
  exposure.ci <- t(apply(retval$exposures, MARGIN = 1, 
                          FUN = quantile, probs = c(prob1, prob2)))
  
  df <- as.data.frame(exposure.ci)
  
  # Get the index where the confidence interval for exposure of a signature contains 0
  idx <- data.table::between(0, lower = df[, 1], upper = df[, 2])
  
  num.of.redundant.sigs <- sum(idx)
  
  # Get rid of redundant sigs and do signature attribution again until we have 
  # no more redundant signatures
  while(num.of.redundant.sigs > 0) {
    # Get the signature names whose confidence interval for exposure does not contain 0
    sig.names <- colnames(sig.universe)[!idx]
    
    # Use sig.names to update the sig.universe and do signature attribution again
    sig.universe <- sig.universe[, sig.names]
    
    retval <- bootstrapSigExposures(m = catalog, 
                                     P = sig.universe, 
                                     R = num.of.replicates,
                                     decomposition.method = method)
    
    exposure.mean <- apply(retval$exposures, MARGIN = 1, FUN = mean)
    
    exposure.ci <- t(apply(retval$exposures, MARGIN = 1, 
                            FUN = quantile, probs = c(prob1, prob2)))
    
    df <- as.data.frame(exposure.ci)
    
    # Get the index where the confidence interval for exposure of a signature contains 0
    idx <- data.table::between(0, lower = df[, 1], upper = df[, 2])
    
    num.of.redundant.sigs <- sum(idx)
  }
  
  exposures.props <- cbind(exposure.mean, exposure.ci)
  exposure.counts <- total.counts * exposures.props
  colnames(exposure.counts)[1] <- colnames(catalog)
  return(exposure.counts)
}

#' Title
#'
#' @inheritParams GetExposureWithConfidence
#' 
#' @param file The name of the PDF file to be produced.
#' 
#' @param ... Optional arguments passed to \code{ICAMSxtra::PlotExposureToPdf}
#'
#' @inheritSection GetExposureWithConfidence Value
#' 
#' @export
GetExposureAndPlotToPdf <- function(catalog, 
                                    file,
                                    sig.universe, 
                                    num.of.replicates = 1000, 
                                    method = decomposeQP,
                                    conf.int = 0.95,
                                    ...) {
  exposure.counts <- GetExposureWithConfidence(catalog = catalog,
                                               sig.universe = sig.universe,
                                               num.of.replicates = num.of.replicates,
                                               method = method,
                                               conf.int = conf.int)
  exposure.counts1 <- exposure.counts[, 1, drop = FALSE]
  ICAMSxtra::PlotExposureToPdf(exposure = exposure.counts1, file = file, ...)
  return(exposure.counts)
}

GetExposureUseSigMiner <- 
  function(catalog, sig.universe, num.of.replicates = 100, conf.int = 0.95) {
    
    retval <- sigminer::sig_fit_bootstrap(catalog = catalog,
                                          sig = sig.universe,
                                          n = num.of.replicates,
                                          method = "QP")
    catalog <- catalog
    total.counts <- colSums(catalog)
    sig.universe <- sig.universe
    num.of.replicates <- num.of.replicates
    conf.int <- conf.int
    exposure.mean <- apply(retval$expo, MARGIN = 1, FUN = mean)
    
    # Use function quantile to find the confidence interval of bootstrap sample values
    prob1 <- (1 - conf.int) / 2
    prob2 <- 1 - prob1
    
    exposure.ci <- t(apply(retval$expo, MARGIN = 1, 
                           FUN = quantile, probs = c(prob1, prob2)))
    
    df <- as.data.frame(exposure.ci)
    
    # Get the index where the confidence interval for exposure of a signature contains 0
    idx <- data.table::between(0, lower = df[, 1], upper = df[, 2])
    
    num.of.redundant.sigs <- sum(idx)
    browser()
    
    # Get rid of redundant sigs and do signature attribution again until we have 
    # no more redundant signatures
    while(num.of.redundant.sigs > 0) {
      # Get the signature names whose confidence interval for exposure does not contain 0
      sig.names <- colnames(sig.universe)[!idx]
      
      # Use sig.names to update the sig.universe and do signature attribution again
      sig.universe <- sig.universe[, sig.names]
      
      sigfit.SBS96.ci <- sigminer::sig_fit_bootstrap(catalog = catalog,
                                                     sig = sig.universe,
                                                     n = num.of.replicates,
                                                     method = "QP")
      
      exposure.mean <- apply(retval$expo, MARGIN = 1, FUN = mean)
      
      exposure.ci <- t(apply(retval$expo, MARGIN = 1, 
                             FUN = quantile, probs = c(prob1, prob2)))
      
      df <- as.data.frame(exposure.ci)
      
      # Get the index where the confidence interval for exposure of a signature contains 0
      idx <- data.table::between(0, lower = df[, 1], upper = df[, 2])
      
      num.of.redundant.sigs <- sum(idx)
    }
    
    exposures.props <- cbind(exposure.mean, exposure.ci)
    sigfit.exposures <- total.counts * exposures.props
    colnames(sigfit.exposures)[1] <- colnames(catalog)
    return(sigfit.exposures)
  }
