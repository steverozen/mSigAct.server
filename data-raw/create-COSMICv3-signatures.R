# Source this file from top level directory

cat(getwd(), "\n")

# Need to have BSgenome.Hsapiens.1000genomes.hs37d5 installed

#################################################################################
# Generate SBS96 signatures for human GRCh37, GRCh37 and GRCm38
file <- "./data-raw/COSMICv3_genome_SBS96.csv"
cat1 <- ICAMS::ReadCatalog(file, ref.genome = "GRCh37",
                         region = "genome",
                         catalog.type = "counts.signature")
COSMIC.v3.hg19.genome.SBS96.sigs <- cat1

COSMIC.v3.hg19.exome.SBS96.sigs <- 
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.SBS96.sigs, 
                          target.ref.genome = "hg19", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.hg38.genome.SBS96.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.SBS96.sigs, 
                          target.ref.genome = "hg38", 
                          target.region = "genome", 
                          target.catalog.type = "counts.signature")
COSMIC.v3.hg38.exome.SBS96.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.SBS96.sigs, 
                          target.ref.genome = "hg38", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.mm10.genome.SBS96.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.SBS96.sigs, 
                          target.ref.genome = "mm10", 
                          target.region = "genome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.mm10.exome.SBS96.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.SBS96.sigs, 
                          target.ref.genome = "mm10", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

##############################################################################
# Generate SBS192 signatures for human GRCh37, GRCh37 and GRCm38
COSMIC.v3.hg19.transcript.SBS192.sigs <- PCAWG7::signature$genome$SBS192

COSMIC.v3.hg19.exome.SBS192.sigs <- 
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.transcript.SBS192.sigs, 
                          target.ref.genome = "hg19", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")
COSMIC.v3.hg38.transcript.SBS192.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.transcript.SBS192.sigs, 
                          target.ref.genome = "hg38", 
                          target.region = "transcript", 
                          target.catalog.type = "counts.signature")
COSMIC.v3.hg38.exome.SBS192.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.transcript.SBS192.sigs, 
                          target.ref.genome = "hg38", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.mm10.transcript.SBS192.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.transcript.SBS192.sigs, 
                          target.ref.genome = "mm10", 
                          target.region = "transcript", 
                          target.catalog.type = "counts.signature")
COSMIC.v3.mm10.exome.SBS192.sigs <-
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.transcript.SBS192.sigs, 
                          target.ref.genome = "mm10", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

##############################################################################
# Generate DBS78 signatures for human GRCh37, GRCh37 and GRCm38

COSMIC.v3.hg19.genome.DBS78.sigs <- PCAWG7::signature$genome$DBS78

COSMIC.v3.hg19.exome.DBS78.sigs <- 
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.DBS78.sigs, 
                          target.ref.genome = "hg19", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.hg38.genome.DBS78.sigs <- 
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.DBS78.sigs, 
                          target.ref.genome = "hg38", 
                          target.region = "genome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.hg38.exome.DBS78.sigs <- 
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.DBS78.sigs, 
                          target.ref.genome = "hg38", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.mm10.genome.DBS78.sigs <- 
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.DBS78.sigs, 
                          target.ref.genome = "mm10", 
                          target.region = "genome", 
                          target.catalog.type = "counts.signature")

COSMIC.v3.mm10.exome.DBS78.sigs <- 
  ICAMS::TransformCatalog(catalog = COSMIC.v3.hg19.genome.DBS78.sigs, 
                          target.ref.genome = "mm10", 
                          target.region = "exome", 
                          target.catalog.type = "counts.signature")

##############################################################################
COSMIC.v3.hg19.genome.ID.sigs <- PCAWG7::signature$genome$ID
ID18.file.path <- "data-raw/sigProfiler_ID_signatures_ID18.csv"
tmp <- ICAMS::ReadCatalog(ID18.file.path, 
                          strict = FALSE, 
                          catalog.type = "counts.signature")
colnames(tmp) <- "ID18"
COSMIC.v3.hg19.genome.ID.sigs <- cbind(COSMIC.v3.hg19.genome.ID.sigs, tmp)


tmp.hg19 <- list(
  genome = list(
    "SBS96"   = COSMIC.v3.hg19.genome.SBS96.sigs,
    "SBS192"   = COSMIC.v3.hg19.transcript.SBS192.sigs,
    "DBS78"  = COSMIC.v3.hg19.genome.DBS78.sigs,
    "ID" = COSMIC.v3.hg19.genome.ID.sigs
  ),
  exome = list(
    "SBS96"   = COSMIC.v3.hg19.exome.SBS96.sigs,
    "SBS192"   = COSMIC.v3.hg19.exome.SBS192.sigs,
    "DBS78"  = COSMIC.v3.hg19.exome.DBS78.sigs
  )
)

tmp.hg38 <- list(
  genome = list(
    "SBS96"   = COSMIC.v3.hg38.genome.SBS96.sigs,
    "SBS192"   = COSMIC.v3.hg38.transcript.SBS192.sigs,
    "DBS78"  = COSMIC.v3.hg38.genome.DBS78.sigs
  ),
  exome = list(
    "SBS96"   = COSMIC.v3.hg38.exome.SBS96.sigs,
    "SBS192"   = COSMIC.v3.hg38.exome.SBS192.sigs,
    "DBS78"  = COSMIC.v3.hg38.exome.DBS78.sigs
  )
)

tmp.mm10 <- list(
  genome = list(
    "SBS96"   = COSMIC.v3.mm10.genome.SBS96.sigs,
    "SBS192"   = COSMIC.v3.mm10.transcript.SBS192.sigs,
    "DBS78"  = COSMIC.v3.mm10.genome.DBS78.sigs
  ),
  exome = list(
    "SBS96"   = COSMIC.v3.mm10.exome.SBS96.sigs,
    "SBS192"   = COSMIC.v3.mm10.exome.SBS192.sigs,
    "DBS78"  = COSMIC.v3.mm10.exome.DBS78.sigs
  )
)

COSMIC.v3.sigs <- list(
  hg19 = tmp.hg19,
  hg38 = tmp.hg38,
  mm10 = tmp.mm10
)