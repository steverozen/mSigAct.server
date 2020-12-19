# Source this file from top level directory

cat(getwd(), "\n")

SBS96.sigs <- COSMIC.v3.genome.SBS96.sigs
for (name in colnames(SBS96.sigs)) {
  file.path <- file.path("inst/app/SBS96", paste0(name, ".png")) 
  grDevices::png(filename=file.path, width = 1700, height = 250)
  ICAMS::PlotCatalog(SBS96.sigs[, name, drop = FALSE])
  grDevices::dev.off()
}

SBS192.sigs <- PCAWG7::signature$genome$SBS192
for (name in colnames(SBS192.sigs)) {
  file.path <- file.path("inst/app/SBS192", paste0(name, ".png")) 
  grDevices::png(filename=file.path, width = 2000, height = 300)
  ICAMS::PlotCatalog(SBS192.sigs[, name, drop = FALSE])
  grDevices::dev.off()
}

DBS78.sigs <- PCAWG7::signature$genome$DBS78
for (name in colnames(DBS78.sigs)) {
  file.path <- file.path("inst/app/DBS78", paste0(name, ".png")) 
  grDevices::png(filename=file.path, width = 2000, height = 300)
  ICAMS::PlotCatalog(DBS78.sigs[, name, drop = FALSE])
  grDevices::dev.off()
}

ID.sigs <- PCAWG7::signature$genome$ID
for (name in colnames(ID.sigs)) {
  file.path <- file.path("inst/app/ID", paste0(name, ".png")) 
  grDevices::png(filename=file.path, width = 1700, height = 250)
  ICAMS::PlotCatalog(ID.sigs[, name, drop = FALSE])
  grDevices::dev.off()
}
