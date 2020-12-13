ui <- fluidPage(
  checkboxGroupInput(inputId = "preselectedSigs",
                     label = HTML('<p>Select metadata <a href=https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/12-2010.1>permacomparisons</a></p>'),
                     choiceNames = foo1,
                     choiceValues = list("SBS2", "SBS3"),
                     selected = c("SBS2", "SBS3")
  ),
  textOutput("txt")
)

server <- function(input, output, session) {
  output$txt <- renderText({
    icons <- paste(input$icons, collapse = ", ")
    paste("You chose", icons)
  })
}

shinyApp(ui, server)


if (FALSE) {
  f1 <- function(){
    shinydashboard::box(title = "SBS1", solidHeader = TRUE,
                        "Propsed aetiology: Spontaneous deamination of 5-methylcytosine (clock-like signature)",
                        img(height = 94, width = 500,src = "www/SBS1.PNG"))
  }
}