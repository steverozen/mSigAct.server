
catSBS192.1 <- PCAWG7::spectra$PCAWG$SBS192[, 1, drop = FALSE]
my.exposure.SBS96 <- PCAWG7::exposure$PCAWG$SBS96[, 1, drop = FALSE]

catDBS78.1 <- PCAWG7::spectra$PCAWG$DBS78[, 1, drop = FALSE]
my.exposure.DBS78 <- PCAWG7::exposure$PCAWG$DBS78[, 1, drop = FALSE]

catID.1 <- PCAWG7::spectra$PCAWG$ID[, 1, drop = FALSE]
my.exposure.ID <- PCAWG7::exposure$PCAWG$ID[, 1, drop = FALSE]


my.sig.SBS192 <-
  CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "SBS192", region = "genome")
my.sig.DBS78 <- 
  CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "DBS78", region = "genome")

my.sig.ID <- 
  CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "ID", region = "genome")


exposures.SBS96 <- GetExposureAndPlotToPdf(catalog = catSBS96.1, 
                                           file = file.path(tempdir(), "test.SBS96.pdf"),
                                           sig.universe = my.sig.SBS96, 
                                           num.of.bootstrap.replicates = 1000, 
                                           method = decomposeQP) 



my.catalog <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
my.sig.universe <- 
  ICAMS.shiny::CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", 
                                     tumor.cohort = "PCAWG",
                                     sig.type = "SBS96", 
                                     region = "genome")
my.exposures <- 
  ICAMS.shiny::GetExposureWithConfidence(catalog = my.catalog,
                                         sig.universe = my.sig.universe,
                                         num.of.bootstrap.replicates = 1000,
                                         method = decomposeQP,
                                         conf.int = 0.95)
my.exposures1 <- 
  ICAMS.shiny::GetExposureAndPlotToPdf(catalog = my.catalog,
                                       file = "D:/test.pdf",
                                       sig.universe = my.sig.universe,
                                       num.of.bootstrap.replicates = 1000,
                                       method = decomposeQP,
                                       conf.int = 0.95)
unlink("Rplots.pdf")

if (FALSE) {
  old.exposure <- findSigExposures(M = my.catalog, 
                                   P = my.sig.universe, 
                                   decomposition.method = decomposeQP)
  old.exposure1 <- bootstrapSigExposures(m = my.catalog,
                                         P = my.sig.universe,
                                         R = 1000, 
                                         decomposition.method = decomposeQP)
}

if (FALSE) {
  sigfit.exposures.SBS96 <- GetExposureUseSigMiner(catalog = catSBS96.1, 
                                                   sig.universe = my.sig.SBS96)
  
  exposures.SBS192 <- GetExposureAndPlotToPdf(catalog = catSBS192.1, 
                                              file = file.path(tempdir(), "test.SBS192.pdf"),
                                              sig.universe = my.sig.SBS192, 
                                              num.of.bootstrap.replicates = 1000, 
                                              method = decomposeQP) 
  sigfit.exposures.SBS192 <- GetExposureUseSigMiner(catalog = catSBS192.1, 
                                                    sig.universe = my.sig.SBS192)
  
  exposures.DBS78 <- GetExposureAndPlotToPdf(catalog = catDBS78.1, 
                                             file = file.path(tempdir(), "test.DBS78.pdf"),
                                             sig.universe = my.sig.DBS78, 
                                             num.of.bootstrap.replicates = 1000, 
                                             method = decomposeQP) 
  
  sigfit.exposures.DBS78 <- GetExposureUseSigMiner(catalog = catDBS78.1, 
                                                   sig.universe = my.sig.DBS78)
  
  exposures.ID <- GetExposureAndPlotToPdf(catalog = catID.1, 
                                          file = file.path(tempdir(), "test.ID.pdf"),
                                          sig.universe = my.sig.ID, 
                                          num.of.bootstrap.replicates = 1000, 
                                          method = decomposeQP) 
  
  sigfit.exposures.ID <- GetExposureUseSigMiner(catalog = catID.1, 
                                                sig.universe = my.sig.ID)
}

