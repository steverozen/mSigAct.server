#' @import shiny
app_server <- function(input, output,session) {
  # List the first level callModules here
  
  # Create an output value that can be used as a condition argument
  # in the function conditionalPanel() in app_ui.R
  output$clicksubmit <- reactive({input$submit})
  outputOptions(output, "clicksubmit", suspendWhenHidden = FALSE) 
  
  # Create an empty list which can be used to store notification ids later
  ids <- list("error" = character(0), "warning" = character(0), 
             "message" = character(0))
  
  # Only when user clicks the submit button, then all the reactive values
  # used as parameters in generating the zip archive will be updated
  observeEvent(
    input$submit,
    {
      output$download <- 
        downloadHandler(filename = paste0(input$zipfile.name, ".zip"),
                        content = function(file) {
                          if (input$vcftype == "strelka.sbs") {
                            # Generate a zip archive from Strelka SBS VCFs and
                            # update the notification ids for errors, warnings
                            # and messages
                            ids <<- 
                              ProcessStrelkaSBSVCFs(input, output, file, ids)
                          } else if (input$vcftype == "strelka.id") {
                            ids <<- ProcessStrelkaIDVCFs(input, output, file, ids)
                          } else if (input$vcftype == "mutect") {
                            ids <<- ProcessMutectVCFs(input, output, file, ids)
                          }
                        })
        
    }
  )
  
  observeEvent(input$submit, {
    #browser()
    RemoveAllNotifications(ids)
  })
}
