#' Transcripts per million (TPM)
#'
#' @rdname tpm
#' @docType methods
#'
#' @param object [bcbioRNADataSet] object.
#'
#' @export
setMethod("tpm", "bcbioRNADataSet", function(object) {
    assays(object)[["tpm"]]
})
