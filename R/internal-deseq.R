#' [DESeqResults] contrast name
#'
#' @rdname res_contrast_name
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
.res_contrast_name <- function(res) {
    mcols(res)[2L, 2L] %>% str_replace("^.*:\\s", "")
}
