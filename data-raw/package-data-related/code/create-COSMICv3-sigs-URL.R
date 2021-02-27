COSMIC.v3.SBS96.sig.names <- colnames(COSMIC.v3.genome.SBS96.sigs)
COSMIC.v3.SBS96.sig.urls <- 
  paste0("https://cancer.sanger.ac.uk/cosmic/signatures/SBS/", 
         COSMIC.v3.SBS96.sig.names, ".tt")
COSMIC.v3.SBS96.sig.links <- 
  matrix(data = COSMIC.v3.SBS96.sig.urls, nrow = length(COSMIC.v3.SBS96.sig.urls), 
         dimnames = list(COSMIC.v3.SBS96.sig.names, "URL"))

COSMIC.v3.SBS192.sig.names <- colnames(COSMIC.v3.genome.SBS96.sigs)
COSMIC.v3.SBS192.sig.urls <- 
  paste0("https://cancer.sanger.ac.uk/cosmic/signatures/SBS/", 
         COSMIC.v3.SBS192.sig.names, ".tt")
COSMIC.v3.SBS192.sig.links <- 
  matrix(data = COSMIC.v3.SBS192.sig.urls, nrow = length(COSMIC.v3.SBS192.sig.urls), 
         dimnames = list(mSigAct:::AddMinusEFor192(COSMIC.v3.SBS96.sig.names), "URL"))

COSMIC.v3.DBS78.sig.names <- paste0("DBS", 1:11)
COSMIC.v3.DBS78.sig.urls <- 
  paste0("https://cancer.sanger.ac.uk/cosmic/signatures/DBS/", 
         COSMIC.v3.DBS78.sig.names, ".tt")
COSMIC.v3.DBS78.sig.links <- 
  matrix(data = COSMIC.v3.DBS78.sig.urls, nrow = length(COSMIC.v3.DBS78.sig.urls), 
         dimnames = list(COSMIC.v3.DBS78.sig.names, "URL"))

COSMIC.v3.ID.sig.names <- paste0("ID", 1:18)
COSMIC.v3.ID.sig.urls <- 
  paste0("https://cancer.sanger.ac.uk/cosmic/signatures/ID/", 
         COSMIC.v3.ID.sig.names, ".tt")
COSMIC.v3.ID.sig.links <- 
  matrix(data = COSMIC.v3.ID.sig.urls, nrow = length(COSMIC.v3.ID.sig.urls), 
         dimnames = list(COSMIC.v3.ID.sig.names, "URL"))
