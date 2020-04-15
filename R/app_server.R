#' @import shiny
app_server <- function(input, output,session) {
  # List the first level callModules here
  
  # Create an empty list which can be used to store notification ids later
  ids <- list("error" = character(0), "warning" = character(0), 
             "message" = character(0))
  
  # Download sample VCFs when user clicks the button
  output$downloadsampleVCFs <- 
    downloadHandler(filename = "sample-VCFs.zip", 
                    content = function(file) {
                      PrepareSampleVCFs(file)
                    })
  
  # Run ICAMS on sample Strelka SBS VCFs when user clicks the button
  output$runstrelkasbsvcfs <- 
    downloadHandler(filename = "ICAMS-test-run-Strelka-SBS-VCFs.zip",
                    content = function(file) {
                      ids <<- 
                        RunICAMSOnSampleStrelkaSBSVCFs(output, file, ids)
                    })
  
  # Run ICAMS on sample Mutect VCFs when user clicks the button
  output$runmutectvcfs <- 
    downloadHandler(filename = "ICAMS-test-run-Mutect-VCFs.zip",
                    content = function(file) {
                      ids <<- 
                        RunICAMSOnSampleMutectVCFs(output, file, ids)
                    })
  
  # When user clicks the submit button, then the program will try to 
  # generate a zip archive based on the input files and parameters
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
        
  # When user clicks the "Remove notifications" button, all the previous
  # notifications(error, warning or message) will be removed
  observeEvent(input$remove, {
    RemoveAllNotifications(ids)
  })
}
