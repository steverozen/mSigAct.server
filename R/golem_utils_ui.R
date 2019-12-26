# Repeat tags$br
#
# @param times the number of br to return
#
# @return the number of br specified in times
# @export
#
# @examples
# rep_br(5)
#
#' @importFrom htmltools HTML
rep_br <- function(times = 1) {
  HTML(rep("<br/>", times = times))
}