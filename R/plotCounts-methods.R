#' @name plotCounts
#' @author Michael Steinbaugh
#' @inherit acidplots::plotCounts
#'
#' @inheritParams acidroxygen::params
#' @param normalized `character(1)`.
#'   Type of normalized counts to visualize.
#' @param ... Additional parameters.
#'
#' @examples
#' data(bcb)
#'
#' g2s <- basejump::Gene2Symbol(bcb)
#' geneIDs <- head(g2s[["geneID"]])
#' print(geneIDs)
#' geneNames <- head(g2s[["geneName"]])
#' print(geneNames)
#'
#' plotCounts(
#'     object = bcb,
#'     genes = geneIDs,
#'     normalized = "vst",
#'     style = "facet"
#' )
#' plotCounts(
#'     object = bcb,
#'     genes = geneNames,
#'     normalized = "vst",
#'     style = "wide"
#' )
NULL



#' @rdname plotCounts
#' @name plotCounts
#' @importFrom bioverbs plotCounts
#' @importMethodsFrom acidplots plotCounts
#' @usage plotCounts(object, ...)
#' @export
NULL



## Updated 2019-07-23.
`plotCounts,bcbioRNASeq` <-  # nolint
    function(object, genes, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)
        ## Ensure counts are always log2 scale.
        if (!normalized %in% c("rlog", "vst")) {
            counts <- log2(counts + 1L)
        }
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts)
        do.call(
            what = plotCounts,
            args = matchArgsToDoCall(
                args = list(
                    object = rse,
                    genes = genes,
                    countsAxisLabel = paste(normalized, "counts (log2)")
                ),
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(`plotCounts,bcbioRNASeq`)
f2 <- methodFormals(
    f = "plotCounts",
    signature = "SummarizedExperiment",
    package = "acidplots"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay", "countsAxisLabel"))]
f <- c(f1, f2)
f[["normalized"]] <- normalizedCounts
formals(`plotCounts,bcbioRNASeq`) <- f



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotCounts",
    signature = signature("bcbioRNASeq"),
    definition = `plotCounts,bcbioRNASeq`
)
