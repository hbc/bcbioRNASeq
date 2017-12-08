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
#' load(system.file(
#'     file.path("inst", "extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' tpm(bcb) %>% head()
NULL



# Methods ====
#' @rdname tpm
#' @export
setMethod(
    "tpm",
    signature("bcbioRNASeq"),
    function(object) {
        assays(object)[["tpm"]]
    })
