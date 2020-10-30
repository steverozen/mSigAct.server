catSBS96.1 <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]

exposures <- GetExposureWithConfidence(catalog = catSBS96.1, 
                                       sig.universe = PCAWG7::signature$genome$SBS96, 
                                       num.of.replicates = 100, 
                                       method = decomposeQP) 

exposure1 <- GetExposureAndPlotToPdf(catalog = catSBS96.1, 
                                     file = file.path(tempdir(), "test.pdf"),
                                     sig.universe = PCAWG7::signature$genome$SBS96, 
                                     num.of.replicates = 100, 
                                     method = decomposeQP) 