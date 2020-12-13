# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

source("data-raw/create-COSMICv3-SBS96-signatures.R")
source("data-raw/create-COSMICv3-sigs-URL.R")
source("data-raw/create-COSMICv3-sigs-aetiology-info.R")

usethis::use_data(COSMIC.v3.genome.SBS96.sigs,
                  COSMIC.v3.SBS.sig.links, 
                  COSMIC.v3.DBS.sig.links,
                  COSMIC.v3.ID.sig.links,
                  SBS.aetiology.HTML,
                  DBS.aetiology.HTML,
                  ID.aetiology.HTML,
                  internal = TRUE, 
                  overwrite = TRUE)
