#' Transcripts Per Million (TPM)
#'
#' @rdname tpm
#' @name tpm
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [matrix].
#'
#' @examples
#' data(bcb)
#' tpm(bcb) %>% head()
NULL



# Methods ====
#' @rdname tpm
#' @export
setMethod(
    "tpm",
    signature("bcbioRNASeqANY"),
    function(object) {
        assays(object)[["tpm"]]
    })
