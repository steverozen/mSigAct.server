#' Run the Shiny Application
#' @param ... Further arguments to be passed to the function.
#'
#' @keywords internal
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(...) {
  with_golem_options(
    app = shinyApp(ui = app_ui, server = app_server, enableBookmarking = "server"),
    golem_opts = list(...)
  )
}

#' @keywords internal
MyRunApp <- function() {
  enableBookmarking(store = "server")
  runApp(list(ui = app_ui, server = app_server),
         host="127.0.0.1", launch.browser = TRUE)
}
