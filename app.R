# Launch the ShinyApp (Do not remove this comment)
# To deploy, run: rsconnect::deployApp()
# Or use the blue button on top of this file
# options("repos" = BiocManager::repositories())
# imports BSgenome.Hsapiens.1000genomes.hs37d5,
# imports BSgenome.Hsapiens.UCSC.hg38,
# imports BSgenome.Mmusculus.UCSC.mm10

# Increase the file uploads limit to 100MB
options(shiny.maxRequestSize = 100*1024^2)
pkgload::load_all(export_all = FALSE,helpers = FALSE,attach_testthat = FALSE)
options( "golem.app.prod" = TRUE)
ICAMS.shiny:::run_app() # add parameters here (if any)