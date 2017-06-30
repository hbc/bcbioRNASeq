#' Sample directories.
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname sample_dirs
#' @docType methods
#'
#' @param object Primary object.
#'
#' @return Folders where samples are kept.
#' @export
setMethod("sample_dirs", "bcbioRNADataSet", function(object) {
    metadata(object)[["sample_dirs"]]
})
