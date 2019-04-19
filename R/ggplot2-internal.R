.dynamicTrans <- function(
    object,
    normalized,
    trans
) {
    validObject(object)
    normalized <- match.arg(normalized, choices = normalizedCounts)
    trans <- match.arg(trans, choices = c("log2", "log10"))

    # Don't allow DESeqTransform to be plotted on log10 scale.
    if (
        normalized %in% c("vst", "rlog") &&
        trans != "log2"
    ) {
        message(paste(normalized, "counts are log2 scale."))
    }

    # Check for DESeqTransform that are already log2 scale and update `trans`.
    if (normalized %in% c("rlog", "vst")) {
        trans <- "identity"
    }

    # Coerce to RangedSummarizedExperiment.
    counts <- counts(object, normalized = normalized)
    rse <- as(object, "RangedSummarizedExperiment")
    assays(rse) <- list(counts)
    assayNames(rse) <- normalized

    # Set the counts axis label.
    countsAxisLabel <- paste(normalized, "counts")

    list(
        object = rse,
        trans = trans,
        countsAxisLabel = countsAxisLabel
    )
}
