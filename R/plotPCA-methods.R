#' @name plotPCA
#' @author Michael Steinbaugh
#' @importMethodsFrom acidplots plotPCA
#' @inherit acidplots::plotPCA
#'
#' @inheritParams params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotPCA(bcb, label = FALSE)
#' plotPCA(bcb, label = TRUE)
NULL



#' @rdname plotPCA
#' @name plotPCA
#' @importFrom BiocGenerics plotPCA
#' @usage plotPCA(object, ...)
#' @export
NULL



## Updated 2019-07-23.
`plotPCA,bcbioRNASeq` <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts."))
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts(object, normalized = normalized))
        do.call(
            what = plotPCA,
            args = matchArgsToDoCall(
                args = list(object = rse),
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(`plotPCA,bcbioRNASeq`)
f2 <- methodFormals(
    f = "plotPCA",
    signature = "SummarizedExperiment",
    package = "acidplots"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
f[["normalized"]] <- normalizedCounts
formals(`plotPCA,bcbioRNASeq`) <- f



#' @rdname plotPCA
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("bcbioRNASeq"),
    definition = `plotPCA,bcbioRNASeq`
)
