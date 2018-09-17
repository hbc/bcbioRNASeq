# FIXME Require DESeqDataSet to use normalized counts



#' Plot Correlation Heatmap
#'
#' This function calculates a correlation matrix based on feature expression per
#' sample.
#'
#' @name plotCorrelationHeatmap
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#' @importFrom basejump plotCorrelationHeatmap
#' @export
#'
#' @inherit basejump::plotCorrelationHeatmap
#'
#' @inheritParams general
#'
#' @seealso
#' - `help("plotCorrelationHeatmap", "basejump")`.
#' - `findMethod("plotCorrelationHeatmap", "SummarizedExperiment")`.
#'
#' @examples
#' plotCorrelationHeatmap(bcb_small, method = "pearson")
#' plotCorrelationHeatmap(bcb_small, method = "spearman")
NULL



.plotCorrelationHeatmap.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle")
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        # Coerce to RangedSummarizedExperiment.
        rse <- as(object, "RangedSummarizedExperiment")
        message(paste("Using", normalized, "counts"))
        counts <- counts(object, normalized = normalized)
        assays(rse) <- list(counts = counts)
        do.call(
            what = plotCorrelationHeatmap,
            args = matchArgsToDoCall(
                args = list(object = rse),
                removeFormals = c("normalized")
            )
        )
    }
f1 <- formals(.plotCorrelationHeatmap.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCorrelationHeatmap",
    signature = "SummarizedExperiment"
)
f <- c(f1, f2[setdiff(names(f2), names(f1))])
formals(.plotCorrelationHeatmap.bcbioRNASeq) <- f



.plotCorrelationHeatmap.DESeqDataSet <-  # nolint
    function(object) {
        validObject(object)
        # Coerce to RangedSummarizedExperiment.
        rse <- as(object, "RangedSummarizedExperiment")
        # Always use normalized counts.
        message("Using normalized counts")
        counts <- counts(object, normalized = TRUE)
        assays(rse) <- list(counts = counts)
        do.call(
            what = plotCorrelationHeatmap,
            args = matchArgsToDoCall(
                args = list(object = rse)
            )
        )
    }
f1 <- formals(.plotCorrelationHeatmap.DESeqDataSet)
f2 <- methodFormals(
    f = "plotCorrelationHeatmap",
    signature = "SummarizedExperiment"
)
f <- c(f1, f2[setdiff(names(f2), names(f1))])
formals(.plotCorrelationHeatmap.DESeqDataSet) <- f



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("bcbioRNASeq"),
    definition = .plotCorrelationHeatmap.bcbioRNASeq
)



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("DESeqDataSet"),
    definition = .plotCorrelationHeatmap.DESeqDataSet
)



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("DESeqTransform"),
    definition = getMethod("plotCorrelationHeatmap", "SummarizedExperiment")
)
