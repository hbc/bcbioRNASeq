#' `DESeqResults` Contrast Name
#'
#' @rdname resContrastName-internal
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
.resContrastName <- function(res) {
    mcols(res)[2L, 2L] %>%
        str_replace("^.*:\\s", "")
}
