# Source this file from top level directory

cat(getwd(), "\n")

file <- "./data-raw/COSMICv3_genome_SBS96.csv"
cat1 <- ICAMS::ReadCatalog(file, ref.genome = "GRCh37",
                         region = "genome",
                         catalog.type = "counts.signature")
original.col.names <- colnames(cat1)
artefacts.sigs.names <- c("SBS27", "SBS43", paste0("SBS", 45:60))
new.col.names <- setdiff(original.col.names, artefacts.sigs.names)

COSMICv3.genome.SBS96.sigs <- cat1[, new.col.names]

usethis::use_data(COSMICv3.genome.SBS96.sigs, internal = TRUE, overwrite = TRUE)
