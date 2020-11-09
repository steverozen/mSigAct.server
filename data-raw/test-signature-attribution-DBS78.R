catDBS78.1 <- PCAWG7::spectra$PCAWG$DBS78[, 1, drop = FALSE]
my.sig.DBS78 <- 
  ICAMS.shiny::CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                                     sig.type = "DBS78", region = "genome")
library(mSigAct)
start.time.DBS78 <- Sys.time()
foo.DBS78 <- SparseAssignActivity(spectra = catDBS78.1, sigs = my.sig.DBS78,
                                  max.level = 8, 
                                  p.thresh = 0.01, 
                                  m.opts = DefaultManyOpts())
end.time.DBS78 <- Sys.time()
time.taken.DBS78 <- end.time.DBS78 - start.time.DBS78
time.taken.DBS78

ground.truth.expo.DBS78 <- PCAWG7::exposure$PCAWG$DBS78[, 1, drop = FALSE]

