  GetExposureUseSigMiner <- 
    function(catalog, sig.universe, num.of.bootstrap.replicates = 100, conf.int = 0.95) {
      
      retval <- sigminer::sig_fit_bootstrap(catalog = catalog,
                                            sig = sig.universe,
                                            n = num.of.bootstrap.replicates,
                                            method = "QP")
      total.counts <- colSums(catalog)
      exposure.mean <- apply(retval$expo, MARGIN = 1, FUN = mean)
      
      # Use function stats::quantile to find the confidence interval of bootstrap sample values
      prob1 <- (1 - conf.int) / 2
      prob2 <- 1 - prob1
      
      exposure.ci <- t(apply(retval$expo, MARGIN = 1, 
                             FUN = stats::quantile, probs = c(prob1, prob2)))
      
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
        
        retval <- sigminer::sig_fit_bootstrap(catalog = catalog,
                                              sig = sig.universe,
                                              n = num.of.bootstrap.replicates,
                                              method = "QP")
        
        exposure.mean <- apply(retval$expo, MARGIN = 1, FUN = mean)
        
        exposure.ci <- t(apply(retval$expo, MARGIN = 1, 
                               FUN = stats::quantile, probs = c(prob1, prob2)))
        
        df <- as.data.frame(exposure.ci)
        
        # Get the index where the confidence interval for exposure of a signature contains 0
        idx <- data.table::between(0, lower = df[, 1], upper = df[, 2])
        
        num.of.redundant.sigs <- sum(idx)
      }
      
      # exposures.props <- cbind(exposure.mean, exposure.ci)
      sigfit.exposures <- cbind(exposure.mean, exposure.ci)
      colnames(sigfit.exposures)[1] <- colnames(catalog)
      return(sigfit.exposures)
    }

