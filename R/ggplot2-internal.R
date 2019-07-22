.dynamicTrans <- function(object, normalized) {
    validObject(object)
    normalized <- match.arg(normalized, choices = normalizedCounts)

    ## Check for DESeqTransform that are already log2 scale and update `trans`.
    if (normalized %in% c("rlog", "vst")) {
        trans <- "identity"
    } else {
        trans <- "log2"
    }

    ## Coerce to RangedSummarizedExperiment.
    counts <- counts(object, normalized = normalized)
    rse <- as(object, "RangedSummarizedExperiment")
    assays(rse) <- list(counts)
    assayNames(rse) <- normalized

    ## Set the counts axis label.
    countsAxisLabel <- paste(normalized, "counts")

    list(
        object = rse,
        trans = trans,
        countsAxisLabel = countsAxisLabel
    )
}
