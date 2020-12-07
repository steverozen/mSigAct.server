# Source this file from top level directory

cat(getwd(), "\n")

SBS96.spectra <- 
  PCAWG7::SplitPCAWGMatrixByTumorType(M = PCAWG7::spectra$PCAWG$SBS96)
SBS96.Lung.AdenoCA.example <- SBS96.spectra$`Lung-AdenoCA`[, 1:2]

SBS192.spectra <- 
  PCAWG7::SplitPCAWGMatrixByTumorType(M = PCAWG7::spectra$PCAWG$SBS192)
SBS192.Lung.AdenoCA.example <- SBS192.spectra$`Lung-AdenoCA`[, 1:2]

SBS1536.spectra <- 
  PCAWG7::SplitPCAWGMatrixByTumorType(M = PCAWG7::spectra$PCAWG$SBS1536)
SBS1536.Lung.AdenoCA.example <- SBS1536.spectra$`Lung-AdenoCA`[, 1:2]

DBS78.spectra <- 
  PCAWG7::SplitPCAWGMatrixByTumorType(M = PCAWG7::spectra$PCAWG$DBS78)
DBS78.Lung.AdenoCA.example <- DBS78.spectra$`Lung-AdenoCA`[, 1:2]

ID.spectra <- 
  PCAWG7::SplitPCAWGMatrixByTumorType(M = PCAWG7::spectra$PCAWG$ID)
ID.Lung.AdenoCA.example <- ID.spectra$`Lung-AdenoCA`[, 1:2]

ICAMS::WriteCatalog(ICAMS::as.catalog(SBS96.Lung.AdenoCA.example), 
                    file = "./inst/extdata/SBS96-mSigAct-example-spectra.csv")

ICAMS::WriteCatalog(ICAMS::as.catalog(SBS192.Lung.AdenoCA.example), 
                    file = "./inst/extdata/SBS192-mSigAct-example-spectra.csv")

ICAMS::WriteCatalog(ICAMS::as.catalog(SBS1536.Lung.AdenoCA.example), 
                    file = "./inst/extdata/SBS1536-mSigAct-example-spectra.csv")

ICAMS::WriteCatalog(ICAMS::as.catalog(DBS78.Lung.AdenoCA.example), 
                    file = "./inst/extdata/DBS78-mSigAct-example-spectra.csv")

ICAMS::WriteCatalog(ICAMS::as.catalog(ID.Lung.AdenoCA.example), 
                    file = "./inst/extdata/ID-mSigAct-example-spectra.csv")
