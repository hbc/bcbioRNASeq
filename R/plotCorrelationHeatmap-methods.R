#' @importFrom basejump plotCorrelationHeatmap
#' @aliases NULL
#' @export
basejump::plotCorrelationHeatmap



#' @name plotCorrelationHeatmap
#' @inherit basejump::plotCorrelationHeatmap
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @examples
#' data(bcb_small)
#' plotCorrelationHeatmap(bcb_small, method = "pearson")
#' plotCorrelationHeatmap(bcb_small, method = "spearman")
NULL



plotCorrelationHeatmap.bcbioRNASeq <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        # Coerce to RangedSummarizedExperiment.
        rse <- as(object, "RangedSummarizedExperiment")
        message(paste("Using", normalized, "counts."))
        counts <- counts(object, normalized = normalized)
        assays(rse) <- list(counts)
        assayNames(rse) <- normalized
        do.call(
            what = plotCorrelationHeatmap,
            args = matchArgsToDoCall(
                args = list(object = rse),
                removeFormals = "normalized"
            )
        )
    }
f1 <- formals(plotCorrelationHeatmap.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCorrelationHeatmap",
    signature = "SummarizedExperiment"
)
f <- c(f1, f2[setdiff(names(f2), c(names(f1), "assay"))])
f[["normalized"]] <- normalizedCounts
formals(plotCorrelationHeatmap.bcbioRNASeq) <- f



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("bcbioRNASeq"),
    definition = plotCorrelationHeatmap.bcbioRNASeq
)
