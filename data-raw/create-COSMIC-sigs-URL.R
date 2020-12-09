
COSMIC.v3.SBS.sig.names <- COSMIC.v3.genome.SBS96.sigs
  #c(colnames(COSMICv3.genome.SBS96.sigs), mSigAct::PossibleArtifacts())
COSMIC.v3.SBS.sig.urls <- 
  paste0("https://cancer.sanger.ac.uk/cosmic/signatures/SBS/", 
         COSMIC.v3.SBS.sig.names, ".tt")
COSMIC.v3.SBS.sig.links <- 
  matrix(data = COSMIC.v3.SBS.sig.urls, nrow = length(COSMIC.v3.SBS.sig.urls), 
         dimnames = list(COSMIC.v3.SBS.sig.names, "URL"))

COSMIC.v3.DBS.sig.names <- paste0("DBS", 1:11)
COSMIC.v3.DBS.sig.urls <- 
  paste0("https://cancer.sanger.ac.uk/cosmic/signatures/DBS/", 
         COSMIC.v3.DBS.sig.names, ".tt")
COSMIC.v3.DBS.sig.links <- 
  matrix(data = COSMIC.v3.DBS.sig.urls, nrow = length(COSMIC.v3.DBS.sig.urls), 
         dimnames = list(COSMIC.v3.DBS.sig.names, "URL"))

COSMIC.v3.ID.sig.names <- paste0("ID", 1:18)
COSMIC.v3.ID.sig.urls <- 
  paste0("https://cancer.sanger.ac.uk/cosmic/signatures/ID/", 
         COSMIC.v3.ID.sig.names, ".tt")
COSMIC.v3.ID.sig.links <- 
  matrix(data = COSMIC.v3.ID.sig.urls, nrow = length(COSMIC.v3.ID.sig.urls), 
         dimnames = list(COSMIC.v3.ID.sig.names, "URL"))
