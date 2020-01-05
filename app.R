# Launch the ShinyApp (Do not remove this comment)
# To deploy, run: rsconnect::deployApp()
# Or use the blue button on top of this file
# options("repos" = BiocManager::repositories())
# imports BSgenome.Hsapiens.1000genomes.hs37d5,
# imports BSgenome.Hsapiens.UCSC.hg38,
# imports BSgenome.Mmusculus.UCSC.mm10

pkgload::load_all()
options( "golem.app.prod" = TRUE)
ICAMS.shiny::run_app() # add parameters here (if any)
