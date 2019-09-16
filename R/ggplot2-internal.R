#' Coerce to SummarizedExperiment with normalized counts as primary assay
#'
#' @note Updated 2019-09-15.
#' @noRd
.normalizedSE <- function(object, normalized) {
    validObject(object)
    normalized <- match.arg(arg = normalized, choices = .normalized)
    message(sprintf("Using %s counts.", normalized))
    counts <- counts(object = object, normalized = normalized)
    assays <- SimpleList(counts)
    names(assays) <- normalized
    se <- as(object, "RangedSummarizedExperiment")
    assays(se) <- assays
    se
}




## FIXME This approach is now broken due to "countsAxisLabel"

#' Dynamic handling of axis transformation
#'
#' Applies to counts that need to be log scaled or are already log2.
#'
#' @note Updated 2019-09-16.
#' @noRd
.dynamicTrans <- function(object, normalized, ...) {
    se <- .normalizedSE(object = object, normalized = normalized)
    dots <- list(...)
    ## FIXME
    assert(areDisjointSets(c("countsAxisLabel", "trans"), names(dots)))
    ## FIXME
    ## > countsAxisLabel <- paste(normalized, "counts")
    ## Check for DESeqTransform that are already log2 scale and update `trans`.
    if (normalized %in% c("rlog", "vst")) {
        trans <- "identity"
        ## FIXME
        ## > countsAxisLabel <- paste(countsAxisLabel, "(log2)")
    } else {
        trans <- "log2"
    }
    ## Return arguments list.
    out <- list(
        object = se,
        trans = trans
        ## FIXME
        ## > countsAxisLabel = countsAxisLabel
    )
    out <- c(out, dots)
    out
}
