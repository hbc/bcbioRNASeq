#' Count matrix accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the `normalized`
#' argument.
#'
#' @rdname counts
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param ... Additional parameters.
#' @param normalized Select raw counts (`FALSE`), DESeq2 normalized counts
#'   (`TRUE`), or additional normalization methods:
#'
#' - `tpm`: Transcripts per million.
#' - `tmm`: Trimmed mean of M-values (edgeR).
#' - `rlog`: Regularized log transformation ([DESeq2::rlog()]).
#' - `vst`: Variance stabilizing transformation
#'   ([DESeq2::varianceStabilizingTransformation()]).
#'
#' @return Counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Raw counts
#' counts(bcb, normalized = FALSE)
#'
#' # DESeq2 normalized counts
#' counts(bcb, normalized = TRUE)
#'
#' # TPM
#' counts(bcb, normalized = "tpm")
#'
#' # TMM
#' counts(bcb, normalized = "tmm")
#'
#' # rlog
#' counts(bcb, normalized = "rlog")
#'
#' # VST
#' counts(bcb, normalized = "vst")
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
        stop("Unsupported normalization method")
    }

    counts <- assays(object)[[slot]]

    # Return matrix from [DESeqTransform]
    if (slot %in% c("rlog", "vst")) {
        counts <- assay(counts)
    }

    counts
})
