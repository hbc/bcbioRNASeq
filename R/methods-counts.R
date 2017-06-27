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
#' or specifically request TMM-normalized counts used for QC functions (`tmm`).
#'
#' @return Counts matrix
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
setMethod("counts", "bcbioRNADataSet", function(object, normalized = FALSE) {
    if (normalized == FALSE) {
        assay(object)
    } else if (isTRUE(normalized) | normalized == "tpm") {
        message("Transcripts per million (TPM)")
        tpm(object)
    } else if  (normalized == "tmm") {
        message("Trimmed mean of M-values normalization")
        tmm(object)
    }  else {
        counts <- assays(object)[[normalized]]
        if (is.null(counts)) {
            stop("Assay missing in SummarizedExperiment")
        }
        message(paste(normalized, "normalization"))
        counts
    }
})
