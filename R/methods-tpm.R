#' Transcripts Per Million (TPM)
#'
#' @rdname tpm
#' @name tpm
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase tpm
#'
#' @inheritParams general
#'
#' @return [matrix].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
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
        data <- assays(object)[["tpm"]]
        assert_is_matrix(data)
        data
    })
