#' Sample directories of [bcbioRnaDataSet] object.
#'
#' This method will be used to access folders where sample information is kept.
#'
#' @rdname sample_dirs
#' @docType methods
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @return Folders where samples are kept.
#' @export
setMethod("sample_dirs", "bcbioRnaDataSet", function(object) {
    metadata(object)[["sample_dirs"]]
})
