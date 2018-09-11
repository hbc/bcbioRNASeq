#' Transcripts Per Million (TPM)
#'
#' @note Both gene- and transcript-level counts are supported, as long as they
#'   were imported using a tximport caller (e.g. salmon, kallisto).
#'
#' @name tpm
#' @family Data Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' x <- tpm(bcb_small)
#' summary(x)
NULL



#' @rdname tpm
#' @export
setMethod(
    "tpm",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)
        assert_is_subset("tpm", assayNames(object))
        assays(object)[["tpm"]]
    }
)
