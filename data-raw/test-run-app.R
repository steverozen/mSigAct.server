require(shiny)
library(DT)

shinyUI(
  DT::dataTableOutput('mytable')
)

# Server.R
library(shiny)
library(DT)


dat <- data.frame(
  country = c('USA', 'China'),
  flag = c('<img src="test.png" height="52"></img>',
           '<img src="http://upload.wikimedia.org/wikipedia/commons/thumb/f/fa/Flag_of_the_People%27s_Republic_of_China.svg/200px-Flag_of_the_People%27s_Republic_of_China.svg.png" height="52"></img>'
  )
)

shinyServer(function(input, output){
  output$mytable <- DT::renderDataTable({
    
    DT::datatable(dat, escape = FALSE) # HERE
  })
})