#' Count Matrix Accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the `normalized`
#' argument.
#'
#' @rdname counts
#' @name counts
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics counts
#'
#' @inheritParams general
#'
#' @param normalized Select raw counts (`FALSE`), DESeq2 normalized counts
#'   (`TRUE`), or additional normalization methods:
#'   - `tpm`: Transcripts per million.
#'   - `tmm`: Trimmed mean of M-values (edgeR).
#'   - `rlog`: Regularized log transformation ([DESeq2::rlog()]).
#'   - `vst`: Variance stabilizing transformation
#'     ([DESeq2::varianceStabilizingTransformation()]).
#'
#' @return [matrix].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # Raw counts (from tximport)
#' counts(bcb, normalized = FALSE) %>% summary()
#'
#' # TPM (from tximport)
#' counts(bcb, normalized = "tpm") %>% summary()
#'
#' # Normalized counts (from DESeq2)
#' counts(bcb, normalized = TRUE) %>% summary()
#'
#' # rlog (from DESeq2)
#' counts(bcb, normalized = "rlog") %>% summary()
#'
#' # TMM (from edgeR)
#' counts(bcb, normalized = "tmm") %>% summary()
#'
#' # VST (from DESeq2)
#' counts(bcb, normalized = "vst") %>% summary()
NULL



# Constructors =================================================================
.counts.bcbioRNASeq <- function(  # nolint
    object,
    normalized = FALSE) {
    if (normalized == FALSE) {
        slot <- "raw"
    } else if (normalized == TRUE) {
        slot <- "normalized"
    } else {
        slot <- normalized
    }

    # Check for slot presence
    if (!slot %in% names(assays(object))) {
         warn(paste(
            slot, "counts not defined.",
            "Calculating and using log2 tmm counts on the fly instead."
        ))
        tmm <- tmm(object)
        return(log2(tmm + 1L))
    }

    counts <- assays(object)[[slot]]

    # Return the matrix from slotted DESeqTransform objects
    if (is(counts, "DESeqTransform")) {
        counts <- assay(counts)
    }

    counts
}


# Methods ======================================================================
#' @rdname counts
#' @export
setMethod(
    "counts",
    signature("bcbioRNASeq"),
    .counts.bcbioRNASeq)
