# Source this file from top level directory

cat(getwd(), "\n")

dir0 <- "inst/app/COSMIC"
dir.create(dir0, showWarnings = FALSE)
for(ref.genome in c("hg19", "hg38", "mm10")) {
  dir1 <- file.path(dir0, ref.genome)
  dir.create(dir1)
  for (region in c("genome", "exome")) {
    dir2 <- file.path(dir1, region)
    dir.create(dir2)
    for (cat.type in c("SBS96", "SBS192", "DBS78")) {
      dir3 <- file.path(dir2, cat.type)
      dir.create(dir3)
      if (cat.type == "SBS96") {
        width <- 1700
        height <- 250
      } else {
        width <- 2000
        height <- 300
      }
      dir <- file.path("inst/app/COSMIC", ref.genome, region, cat.type)
      dir.create(path = dir)
      sigs <- COSMIC.v3.sigs[[ref.genome]][[region]][[cat.type]]
      for (name in colnames(sigs)) {
        file.name <- file.path(dir, paste0(name, ".png"))
        grDevices::png(filename = file.name, width = width, height = height)
        ICAMS::PlotCatalog(sigs[, name, drop = FALSE])
        grDevices::dev.off()
      }
    }
  }
}

ID.sigs <- COSMIC.v3.sigs$hg19$genome$ID
for (name in colnames(ID.sigs)) {
  file.path <- file.path("inst/app/ID", paste0(name, ".png")) 
  grDevices::png(filename=file.path, width = 1700, height = 250)
  ICAMS::PlotCatalog(ID.sigs[, name, drop = FALSE])
  grDevices::dev.off()
}
