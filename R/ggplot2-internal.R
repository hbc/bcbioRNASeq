#' Dynamic handling of axis transformation
#'
#' Applies to counts that need to be log scaled or are already log2.
#'
#' @note Updated 2019-09-15.
#' @noRd
.dynamicTrans <- function(object, normalized, ...) {
    validObject(object)
    normalized <- match.arg(arg = normalized, choices = normalizedCounts)
    dots <- list(...)
    assert(areDisjointSets(c("countsAxisLabel", "trans"), names(dots)))
    countsAxisLabel <- paste(normalized, "counts")
    ## Check for DESeqTransform that are already log2 scale and update `trans`.
    if (normalized %in% c("rlog", "vst")) {
        trans <- "identity"
        countsAxisLabel <- paste(countsAxisLabel, "(log2)")
    } else {
        trans <- "log2"
    }
    ## Coerce to RangedSummarizedExperiment.
    counts <- counts(object, normalized = normalized)
    assays <- SimpleList(counts)
    names(assays) <- normalized
    rse <- as(object, "RangedSummarizedExperiment")
    assays(rse) <- assays
    ## Return arguments list.
    out <- list(
        object = rse,
        trans = trans,
        countsAxisLabel = countsAxisLabel
    )
    out <- c(out, dots)
    out
}
