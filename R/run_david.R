#' DAVID gene set enrichment analysis
#'
#' @author Michael Steinbaugh
#'
#' @importFrom basejump runDavid
#'
#' @param ... Pass through to \code{basejump::runDavid()}
#'
#' @export
run_david <- function(...) {
    basejump::runDavid(...)
}
