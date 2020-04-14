#' @import shiny
app_server <- function(input, output,session) {
  # List the first level callModules here
  
  # Create output values that can be used as condition argument
  # in the function conditionalPanel() in app_ui.R
  output$clicksubmit <- reactive({input$submit})
  output$usebuiltindata <- reactive({input$builtindata})
  outputOptions(output, "clicksubmit", suspendWhenHidden = FALSE) 
  outputOptions(output, "usebuiltindata", suspendWhenHidden = FALSE) 
  
  # Create an empty list which can be used to store notification ids later
  ids <- list("error" = character(0), "warning" = character(0), 
             "message" = character(0))
  
  # Download sample VCFs when user clicks the buttom
  output$downloadsampleVCFs <- downloadHandler(filename = "sample-VCFs.zip", 
                                             content = function(file) {
                                               PrepareSampleVCFs(file)
                                             })
  
  # Generate the zip archive if user chooses to use built-in data
  observeEvent(
    input$builtindata,
    {
      output$download <- 
        downloadHandler(filename = "ICAMS-test.zip",
                        content = function(file) {
                          #browser()
                          ids <<- 
                            GenerateZipFileFromBuiltInData(output, file, ids)
                        })
    })
  
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
  
  # When user clicks the submit buttom, all the previous notifications(error,
  # warning or message) will be removed
  observeEvent(input$submit, {
    RemoveAllNotifications(ids)
  })
}
