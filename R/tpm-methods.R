#' Transcripts Per Million (TPM)
#'
#' @name tpm
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' tpm(bcb_small) %>% summary()
NULL



#' @rdname tpm
#' @export
setMethod(
    "tpm",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)
        assays(object)[["tpm"]]
    }
)
