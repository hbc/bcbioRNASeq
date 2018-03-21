#' Transcripts Per Million (TPM)
#'
#' @name tpm
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase tpm
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' tpm(bcb_small) %>% glimpse()
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
