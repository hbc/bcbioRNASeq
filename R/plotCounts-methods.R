#' @name plotCounts
#' @author Michael Steinbaugh
#' @inherit AcidPlots::plotCounts
#' @note Updated 2019-09-18.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams counts
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#'   See `AcidPlots::plotCounts()` for details.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#'
#' g2s <- basejump::Gene2Symbol(bcb)
#' geneIds <- head(g2s[["geneId"]])
#' print(geneIds)
#' geneNames <- head(g2s[["geneName"]])
#' print(geneNames)
#'
#' plotCounts(bcb, genes = geneIds, style = "facet")
#' plotCounts(bcb, genes = geneNames, style = "wide")
NULL



## Updated 2019-09-18.
`plotCounts,bcbioRNASeq` <-  # nolint
    function(object, genes, normalized, ...) {
        args <- .normalizedPlotArgs(
            object = object,
            genes = genes,
            normalized = match.arg(normalized),
            ...
        )
        do.call(what = plotCounts, args = args)
    }

formals(`plotCounts,bcbioRNASeq`)[["normalized"]] <- .normalized



#' @rdname plotCounts
#' @export
setMethod(
    f = "plotCounts",
    signature = signature("bcbioRNASeq"),
    definition = `plotCounts,bcbioRNASeq`
)
