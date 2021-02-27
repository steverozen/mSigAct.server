# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

source("data-raw/package-data-related/code/create-COSMICv3-signatures.R")
source("data-raw/package-data-related/code/create-COSMICv3-sigs-URL.R")
source("data-raw/package-data-related/code/create-COSMICv3-sigs-aetiology-info.R")

usethis::use_data( COSMIC.v3.hg19.genome.ID.sigs,
                   COSMIC.v3.sigs,
                   COSMIC.v3.SBS96.sig.links, 
                   COSMIC.v3.SBS192.sig.links, 
                   COSMIC.v3.DBS78.sig.links,
                   COSMIC.v3.ID.sig.links,
                   SBS.aetiology,
                   sigs.etiologies,
                   DBS.aetiology,
                   ID.aetiology,
                   internal = TRUE, 
                   overwrite = TRUE)
