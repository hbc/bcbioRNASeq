#' Transcripts per million (TPM)
#'
#' @rdname tpm
#' @docType methods
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @export
setMethod("tpm", "bcbioRnaDataSet", function(object) {
    assays(object)[["abundance"]]
})
