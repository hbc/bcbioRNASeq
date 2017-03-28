#' Import sailfish counts
#'
#' @keywords internal
#' @param ... Passthrough to \code{import_counts()}
#' @export
import_sailfish <- function(...) {
    import_counts(..., type = "sailfish")
}
