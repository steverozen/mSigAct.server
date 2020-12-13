library(shiny)
library(shinyWidgets)

managers <- c('Ram', 'Vijay','Arun','Aswin')
dept <- c('A','B','C','D')
details <- data.frame("Managers" = managers, "Department" = dept, stringsAsFactors = F)

ui <- fluidPage(
  pickerInput(
    'manager', 'Manager',
    choices = managers ,
    c('Ram', 'Vijay','Arun','Aswin'),
    options = list(
      `actions-box` = TRUE),
    multiple = TRUE
  ),
  uiOutput('picker2')
)

server <- function(input, output, session) {
  output$picker2 <- renderUI({
    choices = details$Department[details$Managers %in% input$manager]
    pickerInput('dept', 'Department', choices = choices, choices, options = list(
      `actions-box` = TRUE), multiple = TRUE)
  })
}

shinyApp(ui, server)