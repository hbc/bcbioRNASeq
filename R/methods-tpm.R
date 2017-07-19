#' Transcripts per million (TPM)
#'
#' @rdname tpm
#' @author Michael Steinbaugh
#'
#' @return [matrix].
#' @export
#'
#' @examples
#' data(bcb)
#' tpm(bcb) %>% head
setMethod("tpm", "bcbioRNADataSet", function(object) {
    assays(object)[["tpm"]]
})
