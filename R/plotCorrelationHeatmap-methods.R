#' @name plotCorrelationHeatmap
#' @author Michael Steinbaugh
#' @inherit acidplots::plotCorrelationHeatmap
#' @importMethodsFrom acidplots plotCorrelationHeatmap
#' @note Updated 2019-08-05.
#'
#' @inheritParams plotCounts
#' @inheritParams acidroxygen::params
#'
#' @examples
#' data(bcb)
#' plotCorrelationHeatmap(bcb, method = "pearson")
#' plotCorrelationHeatmap(bcb, method = "spearman")
NULL



#' @rdname plotCorrelationHeatmap
#' @name plotCorrelationHeatmap
#' @importFrom bioverbs plotCorrelationHeatmap
#' @export
NULL



## Updated 2019-07-23.
`plotCorrelationHeatmap,bcbioRNASeq` <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        ## Coerce to RangedSummarizedExperiment.
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

f1 <- formals(`plotCorrelationHeatmap,bcbioRNASeq`)
f2 <- methodFormals(
    f = "plotCorrelationHeatmap",
    signature = "SummarizedExperiment",
    package = "acidplots"
)
f <- c(f1, f2[setdiff(names(f2), c(names(f1), "assay"))])
f[["normalized"]] <- normalizedCounts
formals(`plotCorrelationHeatmap,bcbioRNASeq`) <- f



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("bcbioRNASeq"),
    definition = `plotCorrelationHeatmap,bcbioRNASeq`
)
