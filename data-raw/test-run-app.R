library(shiny)
library(shinydashboard)

#rm(list=ls)

######/ UI Side/######

header <- dashboardHeader(title = "Test")
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("First Tab",tabName = "FTab", icon = icon("globe")),
    menuItem("Second Tab",tabName = "STab", icon = icon("star"))
  ),
  selectInput("navSel", "Selection:", c("b","c"))
)
body <- dashboardBody()

ui <- dashboardPage(header, sidebar, body)


######/ SERVER Side/######

server <- function(input, output, session) {
  
}

shinyApp(ui, server)