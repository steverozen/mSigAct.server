library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(),
  dashboardBody(
    # initialize shinyjs
    shinyjs::useShinyjs(),
    # add custom JS code
    extendShinyjs(text = "shinyjs.hidehead = function(parm){
                                    $('header').css('display', parm);
                                }"),
    actionButton("button","hide header"),
    actionButton("button2","show header")
  )
)

server <- function(input, output) {
  observeEvent(input$button, {
    js$hidehead('none')           
  })
  observeEvent(input$button2, {
    js$hidehead('')           
  })
}

shinyApp(ui, server) 