# Source this file from top level directory

cat(getwd(), "\n")

# Create COSMIC SBS signatures aetiology information
SBS.file <- "data-raw/package-data-related/data/COSMIC-v3-SBS-proposed-aetiology.csv"
SBS.info <- data.table::fread(SBS.file, header = TRUE)
SBS.aetiology <- SBS.info$proposed.aetiology
names(SBS.aetiology) <- SBS.info$name

SBS.aetiology.HTML <-
  paste0('<p><a href="https://cancer.sanger.ac.uk/cosmic/signatures/SBS/',
         rownames(SBS.aetiology), '.tt" target = "_blank">', 
         rownames(SBS.aetiology), '</a> aetiology: ', SBS.aetiology,'</p>')
names(SBS.aetiology.HTML) <- rownames(SBS.aetiology)

# Create COSMIC DBS signatures aetiology information
DBS.file <- "data-raw/package-data-related/data/COSMIC-v3-DBS-proposed-aetiology.csv"
DBS.info <- data.table::fread(DBS.file, header = TRUE)
DBS.aetiology <- DBS.info$proposed.aetiology
names(DBS.aetiology) <- DBS.info$name

DBS.aetiology.HTML <-
  paste0('<p><a href="https://cancer.sanger.ac.uk/cosmic/signatures/SBS/',
         rownames(DBS.aetiology), '.tt" target = "_blank">', 
         rownames(DBS.aetiology), '</a> aetiology: ', DBS.aetiology,'</p>')
names(DBS.aetiology.HTML) <- rownames(DBS.aetiology)

# Create COSMIC ID signatures aetiology information
ID.file <- "data-raw/package-data-related/data/COSMIC-v3-ID-proposed-aetiology.csv"
ID.info <- data.table::fread(ID.file, header = TRUE)
ID.aetiology <- ID.info$proposed.aetiology
names(ID.aetiology) <- ID.info$name

ID.aetiology.HTML <-
  paste0('<p><a href="https://cancer.sanger.ac.uk/cosmic/signatures/SBS/',
         rownames(ID.aetiology), '.tt" target = "_blank">', 
         rownames(ID.aetiology), '</a> aetiology: ', ID.aetiology,'</p>')
names(ID.aetiology.HTML) <- rownames(ID.aetiology)



