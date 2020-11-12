file.ID <- "~/data/PCAWG7-spectrum/BTSG_WGS_PCAWG.indels.csv"
spectrum.ID <- ICAMS::ReadCatalog(file.ID)

file.DBS <- "~/data/PCAWG7-spectrum/BTSG_WGS_PCAWG.dbs.78.csv"
spectrum.DBS <- ICAMS::ReadCatalog(file.DBS)

signatures.ID <- PCAWG7::signature$genome$ID
names.ID <- rownames(foo.ID$exposure)[foo.ID$exposure[,1] >0.00001]

PlotCatalogToPdf(spectrum.ID[, 1, drop = FALSE], file = "~/output/ID.test1.pdf")
PlotCatalogToPdf(spectrum.DBS[, 1, drop = FALSE], file = "~/output/DBS.test1.pdf")


signatures.ID.mSigAct <- 
  
x <- 
  mSigAct:::prop.reconstruct(sigs = signatures.ID[,rownames(foo.ID$exposure)], 
                             exp = foo.ID$exposure)
y <- as.catalog(x)
PlotCatalogToPdf(y, file = "~/output/ID.reconstructed.1.pdf")
PlotCatalogToPdf(cbind(spectrum.ID[, 1, drop = FALSE], zz, y), file = "~/output/compare2.pdf")

z <- mSigAct:::prop.reconstruct(sigs = signatures.ID[,rownames(ground.truth.expo.ID)], 
                               exp = ground.truth.expo.ID)
zz <- as.catalog(z)

foo2 <- philentropy::cosine_dist(spectrum.ID[, 1, drop = FALSE], z, testNA = F)
foo2
foo3 <- philentropy::cosine_dist(spectrum.ID[, 1, drop = FALSE], x, testNA = F)
foo3
