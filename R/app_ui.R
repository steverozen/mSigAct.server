#' @import shiny
#' @import shinythemes
app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here 
    fixedPage(
      # Specify the theme used for UI
      theme = shinytheme("cerulean"),
      
      # Add a title on top the page
      titlePanel(p("ICAMS: In-depth Characterization and Analysis of ",
                 "Mutational Signatures")),
      
      # Add a horizontal line
      hr(),
      
      # Add the first row of control widgets
      fixedRow(
        # Add radio buttons for user to specify the type of VCF files
        column(6, add_vcf_type()),
        
        # Add a conditional panel for user to specify the column names in 
        # Mutect VCFs which contain the tumor sample information (if needed)
        column(6, conditionalPanel(
          condition = "input.vcftype == 'mutect'",
          add_tumor_col_names()))
      ),
      
      # Add the next row of control widgets
      fixedRow(
        # Add radio buttons for user to specify the reference genome
        column(6, add_reference_genome()),
        
        # Add radio buttons for user to specify the genomic region
        # from where the VCFs were generated
        column(6, add_region())
      ),
      
      # Add the next row of control widgets
      fixedRow(
        # Add text input for user to specify the sample names
        # representing different VCF files
        column(6, add_sample_names()),
        
        # Add text input for user to specify the base filename
        # of the CSV and PDF files generated
        column(6, add_base_filename())
      ),
      
      # Add the next row of control widgets
      fixedRow(
        # Add text input for user to specify the zip file name
        column(6, add_zipfile_name()),
        
        # Add a file upload control for user to upload multiple VCF files
        column(6, fileInput(inputId = "vcf.files", label = "Choose VCF Files", 
                            multiple = TRUE))
      ),
      
      # Add one line break
      br(),
      
      # Add an action button for user to update argument and submit
      actionButton(inputId = "submit", label = "Update argument and submit"),
      
      # Add two line breaks
      rep_br(2),
      
      # Add a download button appearing as soon as when user clicks the 
      # action button to update argument and submit for the first time
      conditionalPanel(
        condition = "output.clicksubmit",
        downloadButton(outputId = "download", label = "Download results")),
      
      # Add two line breaks
      rep_br(2),
      
      # Add a horizontal line
      hr(),
      
      # Add a footer on the page showing the URL of ICAMS package
      p("For complete documentation of ICAMS, please refer to ",
        a(href = "https://cran.rstudio.com/web/packages/ICAMS/index.html",
          "https://cran.rstudio.com/web/packages/ICAMS/index.html"))
      
    )
  )
}



#' @import shiny
golem_add_external_resources <- function(){
  
  addResourcePath(
    'www', system.file('app/www', package = 'ICAMS.shiny')
  )
 
  tags$head(
    golem::activate_js(),
    golem::favicon()
    # Add here all the external resources
    # If you have a custom.css in the inst/app/www
    # Or for example, you can add shinyalert::useShinyalert() here
    #tags$link(rel="stylesheet", type="text/css", href="www/custom.css")
  )
}
