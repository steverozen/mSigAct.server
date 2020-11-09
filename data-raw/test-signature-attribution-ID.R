catID.1 <- PCAWG7::spectra$PCAWG$ID[, 1, drop = FALSE]
my.sig.ID <- 
  ICAMS.shiny::CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                                     sig.type = "ID", region = "genome")
library(mSigAct)
start.time.ID <- Sys.time()
foo.ID <- SparseAssignActivity(spectra = catID.1, sigs = my.sig.ID,
                                  max.level = 7, 
                                  p.thresh = 0.01, 
                                  m.opts = DefaultManyOpts())
end.time.ID <- Sys.time()
time.taken.ID <- end.time.ID - start.time.ID
time.taken.ID

ground.truth.expo.ID <- PCAWG7::exposure$PCAWG$ID[, 1, drop = FALSE]
