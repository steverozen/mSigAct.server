library(shiny)
mymtcars = mtcars
mymtcars$id = 1:nrow(mtcars)
runApp(
  list(ui = pageWithSidebar(
    headerPanel('Examples of DataTables'),
    sidebarPanel(
      shinyWidgets::pickerInput(inputId = 'selectedSig', label = 'Selected signatures', 
                         choices = c("SBS1", "SBS2"),
                         selected = c("SBS1", "SBS2"),
                         options = list(
                           `actions-box` = TRUE), 
                         multiple = TRUE),
      textInput("collection_txt",label="Foo")
    ),
    mainPanel(
      dataTableOutput("mytable")
    )
  )
  , server = function(input, output, session) {
    
    dat <- data.frame(
      name = c('<a href="http://rstudio.com">SBS1</a>', '<a href="http://rstudio.com">SBS2</a>'),
      spectrum = c('<img src="SBS/SBS1.png" height="52"></img>',
               '<img src="SBS/SBS2.png" height="52"></img>'),
      proposed.aetiology = c('Spontaneous deamination of 5-methylcytosine (clock-like signature)',
                    'Activity of APOBEC family of cytidine deaminases')
    )
    
    rownames(dat) <- c("SBS1", "SBS2")
    
    
    output$mytable <- DT::renderDataTable({
      
      DT::datatable(dat[input$selectedSig, ], escape = FALSE, rownames = FALSE) # HERE)
  })
  })
)