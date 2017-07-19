#' `DESeqResults` Contrast Name
#'
#' @rdname res_contrast_name
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @param res [DESeqResults].
#'
#' @return Contrast name string.
.res_contrast_name <- function(res) {
    mcols(res)[2L, 2L] %>% str_replace("^.*:\\s", "")
}
