#' @importFrom purrr set_names
#' @keywords internal
#' @export
set_names <- function(...) {
    purrr::set_names(...)
}



#' import basejump
#' @keywords internal
#' @export
set_rownames <- function(...) {
    basejump::setRownames(...)
}
