ui <- fluidPage(
  checkboxGroupInput(inputId = "preselectedSigs",
                                    label = paste0("These signatures were preselected based ",  
                                                   "on cancer type."),
                                    choiceNames = 
                                      list(HTML('<img height="94" width="500" src="www/SBS1.PNG"/>'),
                                           HTML('<img height="94" width="500" src="www/SBS2.PNG"/>'),
                                           shinydashboard::box("SBS1",
                                                               img(height = 94, width = 500,src = "www/SBS1.PNG")),
                                           shinydashboard::box(
                                             title = "Inputs", solidHeader = TRUE,
                                             "Box content here", br(), "More box content",
                                             sliderInput("slider", "Slider input:", 1, 100, 50),
                                             textInput("text", "Text input:")
                                           )),
                                    choiceValues = list("SBS1", "SBS2", "SBS3", "SBS4"),
                                    selected = c("SBS1", "SBS2", "SBS3", "SBS4")
                                    
  )
 
  ,
  textOutput("txt")
)

server <- function(input, output, session) {
  output$txt <- renderText({
    icons <- paste(input$icons, collapse = ", ")
    paste("You chose", icons)
  })
}

shinyApp(ui, server)

