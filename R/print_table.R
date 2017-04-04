#' Print table
#'
#' Wrapper function for knit reports
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom basejump printTable
#'
#' @param ... Pass through to \code{knitr::kable()}
#'
#' @export
print_table <- function(...) {
    basejump::printTable(...)
}
