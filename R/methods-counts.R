#' Count matrix accessors
#'
#' @rdname counts
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request TMM-normalized counts used for QC functions (`tmm`).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Raw counts
#' counts(bcb, normalized = FALSE)
#'
#' # TPM
#' counts(bcb, normalized = TRUE)
#' counts(bcb, normalized = "tpm")
#'
#' # TMM
#' counts(bcb, normalized = "tmm")
#' }
setMethod("counts", "bcbioRnaDataSet", function(object, normalized = FALSE) {
    if (normalized == "tmm") {
        message("TMM-normalized counts")
        tmm(object)
    } else if (isTRUE(normalized) | normalized == "tpm") {
        message("Transcripts per million (TPM)")
        tpm(object)
    } else {
        message("Raw counts")
        assay(object)
    }
})
