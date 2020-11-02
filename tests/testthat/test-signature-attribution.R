catSBS96.1 <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
my.exposure.SBS96 <- PCAWG7::exposure$PCAWG$SBS96[, 1, drop = FALSE]
catSBS192.1 <- PCAWG7::spectra$PCAWG$SBS192[, 1, drop = FALSE]

my.sig.SBS96 <- 
  CancerTypeToSigSubset(ca.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "SBS96", region = "genome")
my.sig.SBS192 <-
  CancerTypeToSigSubset(ca.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "SBS192", region = "genome")

exposures.SBS96 <- GetExposureAndPlotToPdf(catalog = catSBS96.1, 
                                     file = file.path(tempdir(), "test.SBS96.pdf"),
                                     sig.universe = my.sig.SBS96, 
                                     num.of.replicates = 1000, 
                                     method = decomposeQP) 

exposures.SBS192 <- GetExposureAndPlotToPdf(catalog = catSBS192.1, 
                                            file = file.path(tempdir(), "test.SBS192.pdf"),
                                            sig.universe = my.sig.SBS192, 
                                            num.of.replicates = 1000, 
                                            method = decomposeQP) 


















###
sigfit.SBS96 <- sigminer::sig_fit(catalogue_matrix = catSBS96.1,
                                  sig = my.sig.SBS96,
                                  method = "QP")
sigfit.SBS96.ci <- sigminer::sig_fit_bootstrap(catalog = catSBS96.1,
                                               sig = my.sig.SBS96,
                                               n = 1000,
                                               method = "QP")

sigfit.exposures <- GetExposureUseSigMiner(catalog = catSBS96.1, 
                                           sig.universe = my.sig.SBS96,
                                           num.of.replicates = 100,
                                           conf.int = 0.95)
