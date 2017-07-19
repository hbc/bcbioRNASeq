#' Sample Directories
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname sample_dirs
#' @author Michael Steinbaugh
#'
#' @return Named character vector containing sample directory paths.
#' @export
#'
#' @examples
#' data(bcb)
#' sample_dirs(bcb)
setMethod("sample_dirs", "bcbioRNADataSet", function(object) {
    metadata(object)[["sample_dirs"]]
})
