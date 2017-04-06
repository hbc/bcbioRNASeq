#' Knit table
#'
#' Wrapper function for knit reports
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @importFrom basejump knitTable
#'
#' @param ... Pass through to \code{basejump::knitTable()}
#'
#' @export
knit_table <- function(...) {
    basejump::knitTable(...)
}
