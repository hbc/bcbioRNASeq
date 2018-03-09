#' Transcripts Per Million (TPM)
#'
#' @name tpm
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase tpm
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' tpm(bcb) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname tpm
#' @export
setMethod(
    "tpm",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)
        data <- assays(object)[["tpm"]]
        assert_is_matrix(data)
        data
    }
)
