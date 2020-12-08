# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

source("data-raw/create-COSMICv3-SBS96-signatures.R")
source("data-raw/create-COSMIC-sigs-URL.R")

usethis::use_data(COSMIC.v3.genome.SBS96.sigs,
                  COSMIC.v3.SBS.sig.links, 
                  COSMIC.v3.DBS.sig.links,
                  COSMIC.v3.ID.sig.links,
                  internal = TRUE, 
                  overwrite = TRUE)
