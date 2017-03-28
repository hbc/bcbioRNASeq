#' Import salmon counts
#'
#' @keywords internal
#' @param ... Passthrough to \code{import_counts()}
#' @export
import_salmon <- function(...) {
    import_counts(..., type = "salmon")
}
