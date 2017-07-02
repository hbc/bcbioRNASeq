#' Count matrix accessors
#'
#' By default, [counts()] returns the raw counts. This method will return
#' transcripts per million (TPM) by default when `normalized = TRUE`. This can
#' be overriden by requesting the normalization method directly.
#'
#' @rdname counts
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param ... Additional parameters.
#' @param normalized Select normalized counts (`TRUE`), raw counts (`FALSE`),
#' or specifically request additional normalization methods:
#'
#' - `tpm`: Transcripts per million.
#' - `tmm`: Trimmed mean of M-values.
#'
#' @return Counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Raw counts
#' counts(bcb, normalized = FALSE)
#'
#' # Normalized counts
#' counts(bcb, normalized = TRUE)
#'
#' # TPM
#' counts(bcb, normalized = "tpm")
#'
#' # TMM
#' counts(bcb, normalized = "tmm")
#' }
setMethod("counts", "bcbioRNADataSet", function(object, normalized = FALSE) {
    if (normalized == FALSE) {
        slot <- "raw_counts"
    } else if (normalized == TRUE) {
        slot <- "normalized_counts"
    } else {
        slot <- normalized
    }

    # Check for slot presence
    if (!slot %in% names(assays(object))) {
    }

    assays(object)[[slot]]
})
