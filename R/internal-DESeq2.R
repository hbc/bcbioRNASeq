#' `DESeqResults` Contrast Name
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom S4Vectors mcols
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
.resContrastName <- function(res) {
    mcols(res)[2, 2] %>%
        gsub(x = ., pattern = "^.*:\\s", replacement = "") %>%
        gsub(x = ., pattern = "_", replacement = " ") %>%
        # Improve appearance for difference of differences
        gsub(x = ., pattern = "\\+", replacement = " \\+\n    ")
}
