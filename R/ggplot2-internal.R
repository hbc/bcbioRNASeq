## Dynamic handling of axis transformation.
## Applies to counts that need to be log scaled or are already log2.
## Updated 2019-08-20.
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
    ## Return.
    list(
        object = rse,
        trans = trans,
        countsAxisLabel = countsAxisLabel
    )
}
