#' `DESeqResults` Contrast Name
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
.resContrastName <- function(res) {
    mcols(res)[2, 2] %>%
        str_replace("^.*:\\s", "") %>%
        str_replace_all("_", " ") %>%
        # Improve appearance for difference of differences
        str_replace_all("\\+", " \\+\n    ")
}
