###############################################################
UI
##############################################################

if (FALSE) {
  navlistPanel(id = "navlistResults",
               tabPanel(title = "Best result",
                        tabsetPanel(
                          tabPanel(title = "Attribution counts", 
                                   value = "attributionCountsBest",
                                   uiOutput(outputId = "exposureTable")),
                          tabPanel(title = "Attribution plot", 
                                   value = "attributionPlotBest",
                                   uiOutput(outputId = "pdfview"))),
                        #tabPanel(title = "Signature presence test",
                        #         value = "sigPresenceTestBest")
                        value = "bestResult"),
               tabPanel(title = "Second best result",
                        tabsetPanel(
                          tabPanel(title = "Attribution counts", 
                                   value = "attributionCountsSecond",
                                   uiOutput(outputId = "exposureTable2")),
                          tabPanel(title = "Attribution plot", 
                                   value = "attributionPlotSecond",
                                   uiOutput(outputId = "pdfview2"))),
                        #tabPanel(title = "Signature presence test",
                        #         value = "sigPresenceTestSecond")
                        value = "secondBestResult"),
               tabPanel(title = "Third best result",
                        tabsetPanel(
                          tabPanel(title = "Attribution counts", 
                                   value = "attributionCountsThird",
                                   uiOutput(outputId = "exposureTable3")),
                          tabPanel(title = "Attribution plot", 
                                   value = "attributionPlotThird",
                                   uiOutput(outputId = "pdfview3"))),
                        #tabPanel(title = "Signature presence test",
                        #         value = "sigPresenceTestThird")
                        value = "thirdBestResult"),
               widths = c(2, 10), fluid = FALSE)
}

###############################################################
Server
##############################################################
#' @importFrom dplyr bind_rows
#' @keywords internal
PrepareAttributionResults2 <- 
  function (input, output, session, input.catalog.type, plotdata) {
    
    spect <- plotdata$spect
    sig.universe <- plotdata$sig.universe
    #cossim <- plotdata$cossim
    retval <- plotdata$retval
    
    all.tested <- retval$all.tested
    
    # Sort the rows according to MAP value
    all.tested <- dplyr::arrange(all.tested, desc(MAP))
    
    results.to.show <- 3
    available.results <- min(results.to.show, nrow(all.tested))
    results <- list()
    
    
    for (i in 1:available.results) {
      local({
        my_i <- i
        tabpanelname <- paste0("tabpanel", my_i)
        exposureTablename <- paste0("exposureTable", my_i)
        pdfviewname <- paste0("pdfview", my_i)
        tabsetpanelname <- paste0("result", my_i)
        
        
        result <- all.tested[i, ]
        exp <- result$exp[[1]]
        MAP.exp <- dplyr::tibble(sig.id = names(exp), count = exp)
        
        QP.exp <- 
          mSigAct::OptimizeExposureQP(spect, 
                                      sig.universe[ , MAP.exp$sig.id, 
                                                    drop = FALSE])
        QP.MAP.exp <-
          dplyr::tibble(sig.id = names(QP.exp), count = QP.exp)
        
        r.qp <- mSigAct::ReconstructSpectrum(sig.universe, exp = QP.exp, 
                                             use.sig.names = TRUE)
        reconstructed.catalog0 <- 
          as.catalog(r.qp, ref.genome = input$ref.genome2, 
                     region = input$region2)
        
        cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
        
        colnames(reconstructed.catalog0) <- 
          paste0("reconstructed (cosine similarity = ", cossim, ")")
        reconstructed.catalog <- round(reconstructed.catalog0)
        
        
        tbl1 <- data.frame(names = colnames(spect), count = colSums(spect), 
                           cosine.similarity = cossim)
        tbl2 <- data.frame(names = QP.MAP.exp$sig.id, 
                           count = QP.MAP.exp$count)
        tbl <- dplyr::bind_rows(tbl1, tbl2)
        
        output[[exposureTablename]] <- renderTable({
          tbl
        }, sanitize.text.function = function(x) x, digits = 5)
        
        output[[tabpanelname]] <- renderUI(
          {
            tabPanel(
              title = "Best result",
              tabsetPanel(
                tabPanel(title = "Attribution counts", 
                         #value = "attributionCounts",
                         uiOutput(outputId = exposureTablename)),
                tabPanel(title = "Attribution plot", 
                         #value = "attributionPlot",
                         uiOutput(outputId = pdfviewname)))
              #value = tabsetpanelname)
              
            )
          }
        )
        
      })
    }
    
    
    tabpanel.list <- lapply(1:available.results, function(i) {
      tabpanelname <- paste0("tabpanel", i)
      
      uiOutput(outputId = tabpanelname)
      
    })
    
    do.call(tagList, tabpanel.list)
    
    
    # TODO: The navlist Panel is not working properly, need to find a way
    output$attributionResults2 <- renderUI({
      tagList(
        navlistPanel(id = "navlistResults",
                     tabpanel.list,
                     widths = c(2, 10), fluid = FALSE)
      )
    })
    
    # Show the new attribution results
    shinyjs::show(selector = '#panels li a[data-value=attributionResultsTab]')
    
    shinydashboard::updateTabItems(session = session, inputId = "panels", 
                                   selected = "attributionResultsTab")
  }