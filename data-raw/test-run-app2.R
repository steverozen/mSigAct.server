library(DT)

ui <- basicPage(
  h2("The mtcars data"),
  DT::dataTableOutput("mytable")
)

server <- function(input, output) {
  dat <- data.frame(
    country = c('USA', 'China'),
    flag = c('<img src="DBS78/DBS1.PNG" height="52"></img>',
             '<img src="http://upload.wikimedia.org/wikipedia/commons/thumb/f/fa/Flag_of_the_People%27s_Republic_of_China.svg/200px-Flag_of_the_People%27s_Republic_of_China.svg.png" height="52"></img>'
    )
  )
  output$mytable <- DT::renderDataTable({
    
    DT::datatable(dat, escape = FALSE) # HERE
  })
}

shinyApp(ui, server)