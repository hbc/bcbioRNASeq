#' Transcripts per million (TPM)
#'
#' @rdname tpm
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#'
#' @export
setMethod("tpm", "bcbioRNADataSet", function(object) {
    assays(object)[["tpm"]]
})
