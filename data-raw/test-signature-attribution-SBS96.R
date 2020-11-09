catSBS96.1 <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
my.sig.SBS96 <- 
  ICAMS.shiny::CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                        sig.type = "SBS96", region = "genome")
library(mSigAct)

foo <- SparseAssignActivity(spectra = catSBS96.1, sigs = my.sig.SBS96,
                            max.level = 9, 
                            p.thresh = 0.01, 
                            m.opts = DefaultManyOpts())
start.time.SBS96 <- Sys.time()
end.time.SBS96 <- Sys.time()
time.taken.SBS96 <- end.time.SBS96 - start.time.SBS96
time.taken.SBS96

my.sig.SBS96.2 <- my.sig.SBS96[, 1:17, drop = FALSE]
foo2 <- SparseAssignActivity(spectra = catSBS96.1, sigs = my.sig.SBS96.2,
                            max.level = 9, 
                            p.thresh = 0.01, 
                            m.opts = DefaultManyOpts()) 

my.sig.SBS96.3 <- my.sig.SBS96[, 1:16, drop = FALSE]
foo3 <- SparseAssignActivity(spectra = catSBS96.1, sigs = my.sig.SBS96.3,
                             max.level = 16, 
                             p.thresh = 0.01, 
                             m.opts = DefaultManyOpts()) 

my.sig.SBS96.4 <- my.sig.SBS96[, 1:15, drop = FALSE]
foo4 <- SparseAssignActivity(spectra = catSBS96.1, sigs = my.sig.SBS96.4,
                             max.level = 15, 
                             p.thresh = 0.01, 
                             m.opts = DefaultManyOpts()) 
