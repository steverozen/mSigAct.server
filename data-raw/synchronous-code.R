submitAttribution <- reactive({
  list(input$submitAttributionOnTop, input$submitAttribution2)
})

observeEvent(
  submitAttribution(),
  {
    if(is.null(input$submitAttributionOnTop) && is.null(input$submitAttribution2)){
      return()
    }
    
    if(input$submitAttributionOnTop == 0 && input$submitAttribution2 == 0){
      return()
    }
    
    #Don't do anything if in the middle of a run
    if(running())
      return(NULL)
    running(TRUE)
    
    
    if (length(plot.names) > 0) {
      for (i in 1:length(plot.names)) {
        shinyjs::hide(id = plot.names[i])
      }
    }
    
    spect <- catalog[, input$selectedSampleFromCatalogForAttribution, drop = FALSE]
    catalog.type <- input$selectedCatalogType
    cancer.type <- input$selectedCancerType
    region <- input$region2
    
    if (catalog.type == "SBS192") {
      sig.universe <- 
        PCAWG7::signature[["genome"]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- 
        PCAWG7::exposure.stats$PCAWG[["SBS96"]][[cancer.type]][colnames(sig.universe), ]
    } else {
      sig.universe <- 
        PCAWG7::signature[[region]][[catalog.type]][, input$selectedSigSubset2]
      sigs.prop <- 
        PCAWG7::exposure.stats$PCAWG[[catalog.type]][[cancer.type]][colnames(sig.universe), ]
    }
    
    sig.names <- rownames(sigs.prop)
    sigs.prop <- unlist(sigs.prop[ , 2])
    names(sigs.prop) <- sig.names
    
    if (FALSE) {
      QP.exposure <- 
        GetExposureWithConfidence(catalog = spect,
                                  sig.universe = sig.universe,
                                  num.of.bootstrap.replicates = 10000,
                                  method = decomposeQP,
                                  conf.int = 0.95)
      updated.sig.universe <- sig.universe[ , rownames(QP.exposure)]
      updated.sigs.prop <- sigs.prop[rownames(QP.exposure)]
    }
    
    # Create a Progress object
    progress <- ipc::AsyncProgress$new(session, min = 0, max = 1,
                                       message = "Analysis in progress",
                                       detail = "This may take a while...")
    result_val(NULL)
    
    fut <- future::future(
      {
        # Close the progress when this reactive exits (even if there's an error)
        # on.exit(progress$close())
        
        # Create a callback function to update progress. Each time this is called, it
        # will increase the progress by that value and update the detail
        updateProgress <- function(value = NULL, detail = NULL) {
          
          # TODO: Need to change the callback function in mSigAct::MAPAssignActivityInternal
          # each.level.callback.fn(
          # value = 1/max.level - 0.01,
          # detail = paste0("Testing removal of subsets of ", df, " signatures (",
          #                 length(subsets2), " subsets)"))
          progress$inc(amount = value, detail = detail)
          interruptor$execInterrupts()
        }
        
        retval <- mSigAct::MAPAssignActivity1(
          spect = spect,
          sigs = sig.universe,
          sigs.presence.prop = sigs.prop,
          max.level = length(sigs.prop) - 1,
          p.thresh = 0.01,
          eval_f = mSigAct::ObjFnBinomMaxLHRound,
          eval_g_ineq = mSigAct::g_ineq_for_ObjFnBinomMaxLH2,
          m.opts = mSigAct::DefaultManyOpts(),
          max.mc.cores = 50,
          progress.monitor = updateProgress
        )
        
        return(retval)
      }, seed = TRUE) %...>% {
        retval <- .
        
        
        if (retval$success == FALSE || is.null(retval$success)) {
          output$sigContributionPlot <- renderUI({
            output$attributionMessage <- 
              renderText(paste0("The algorithm could not find the optimal number of ", 
                                "signatures that explain the spectrum. Please reduce the ", 
                                "number of signatures used."))
            tagList(
              textOutput(outputId = "attributionMessage")
            )
          })
        } else {
          MAP.best.exp <- retval$MAP
          
          QP.exp <- 
            mSigAct::OptimizeExposureQP(spect, 
                                        sig.universe[ , MAP.best.exp$sig.id, 
                                                      drop = FALSE])
          QP.best.MAP.exp <-
            tibble::tibble(sig.id = names(QP.exp), QP.best.MAP.exp = QP.exp)
          
          r.qp <- mSigAct::ReconstructSpectrum(sig.universe, exp = QP.exp, 
                                               use.sig.names = TRUE)
          reconstructed.catalog0 <- 
            as.catalog(r.qp, ref.genome = input$ref.genome2, 
                       region = input$region2)
          
          cossim <- round(mSigAct::cossim(spect, reconstructed.catalog0), 5)
          
          colnames(reconstructed.catalog0) <- 
            paste0("reconstructed (cosine similarity = ", cossim, ")")
          reconstructed.catalog <- round(reconstructed.catalog0)
          
          plotdata$spect <<- spect
          plotdata$reconstructed.catalog <<- reconstructed.catalog
          plotdata$sig.universe <<- sig.universe
          plotdata$QP.best.MAP.exp <<- QP.best.MAP.exp
          
          max_plots <- nrow(QP.best.MAP.exp) + 2
          output$sigContributionPlot <- renderUI({
            plot_output_list <- lapply(1:max_plots, function(i) {
              plotname <- paste("plot", i, sep="")
              plot.names[i] <<- plotname
              plotOutput(plotname)
            })
            
            tagList(plot_output_list)
          })
          
          for (i in 1:max_plots) {
            # Need local so that each item gets its own number. Without it, the value
            # of i in the renderPlot() will be the same across all instances, because
            # of when the expression is evaluated.
            local({
              my_i <- i
              plotname <- paste("plot", my_i, sep="")
              
              if (my_i == 1) {
                output[[plotname]] <- renderPlot(
                  expr = ICAMS::PlotCatalog(spect)
                  , width = 800, height = 260
                )
              } else if (my_i == 2) {
                output[[plotname]] <- renderPlot(
                  expr = ICAMS::PlotCatalog(reconstructed.catalog)
                  , width = 800, height = 260
                )
              } else {
                output[[plotname]] <- renderPlot({
                  sig.name <- QP.best.MAP.exp$sig.id[my_i-2]
                  sig.catalog <- sig.universe[, sig.name, drop = FALSE]
                  colnames(sig.catalog) <- 
                    paste0(sig.name, " (exposure = ", 
                           round(QP.best.MAP.exp$QP.best.MAP.exp[my_i-2]), ")")
                  ICAMS::PlotCatalog(sig.catalog)
                  
                }, width = 800, height = 260)
              }
            })
          }
          
          for (i in length(plot.names)) {
            shinyjs::show(id = plot.names[i])
          }
          
        }
        
      } %...>% result_val
    
    # Show notification on error or user interrupt
    fut <- promises::catch(fut,
                           function(e){
                             result_val(NULL)
                             print(e$message)
                             showNotification(e$message)
                           })
    
    # When done with analysis, remove progress bar
    fut <- promises::finally(fut, function(){
      progress$close()
      running(FALSE) # Declare done with run
    })
    
    # Return something other than the future so we don't block the UI
    NULL
  })