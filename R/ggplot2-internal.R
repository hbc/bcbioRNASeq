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



#' Coerce to SummarizedExperiment and dynamically handle axis transformation
#'
#' Applies to counts that need to be log scaled or are already log2.
#'
#' @note Updated 2019-09-16.
#' @noRd
.normalizedPlotArgs <- function(object, normalized, ...) {
    se <- .normalizedSE(object = object, normalized = normalized)
    dots <- list(...)
    assert(areDisjointSets("trans", names(dots)))
    args <- c(
        object = se,
        trans = "log2",
        dots
    )
    countAxis <- paste(normalized, "counts")
    ## Check for DESeqTransform that are already log2 scale and update `trans`.
    if (normalized %in% c("rlog", "vst")) {
        args[["trans"]] <- "identity"
        countAxis <- paste(countAxis, "(log2)")
    }
    ## Automatically stash counts axis label.
    labels <- args[["labels"]]
    if (is.null(labels)) {
        labels <- list()
    }
    labels[["countAxis"]] <- countAxis
    args[["labels"]] <- labels
    args
}
