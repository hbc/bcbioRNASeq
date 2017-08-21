#' `DESeqResults` Contrast Name
#'
#' @rdname internal-resContrastName
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
.resContrastName <- function(res) {
    mcols(res)[2L, 2L] %>%
        str_replace("^.*:\\s", "") %>%
        str_replace_all("_", " ") %>%
        # Improve appearance for difference of differences
        str_replace_all("\\+", " \\+\n    ")
}
