#' @importFrom basejump plotGene
#' @aliases NULL
#' @export
basejump::plotGene



#' @name plotGene
#' @inherit basejump::plotGene
#' @author Michael Steinbaugh
#'
#' @inheritParams basejump.globals::params
#'
#' @examples
#' data(bcb_small, deseq_small)
#'
#' ## bcbioRNASeq ====
#' object <- bcb_small
#' g2s <- Gene2Symbol(object)
#' geneIDs <- head(g2s[["geneID"]])
#' print(geneIDs)
#' geneNames <- head(g2s[["geneName"]])
#' print(geneNames)
#'
#' plotGene(
#'     object = object,
#'     genes = geneIDs,
#'     normalized = "vst",
#'     style = "facet"
#' )
#' plotGene(
#'     object = object,
#'     genes = geneNames,
#'     normalized = "vst",
#'     style = "wide"
#' )
#'
#' ## DESeqAnalysis ====
#' object <- deseq_small
#' plotGene(object, genes = geneIDs, style = "facet")
#' plotGene(object, genes = geneNames, style = "wide")
NULL



plotGene.bcbioRNASeq <-  # nolint
    function(
        object,
        genes,
        normalized
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)
        # Ensure counts are always log2 scale.
        if (!normalized %in% c("rlog", "vst")) {
            counts <- log2(counts + 1L)
        }
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts)
        do.call(
            what = plotGene,
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
f1 <- formals(plotGene.bcbioRNASeq)
f2 <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f2 <- f2[setdiff(names(f2), c(names(f1), "assay", "countsAxisLabel"))]
f <- c(f1, f2)
f[["normalized"]] <- normalizedCounts
formals(plotGene.bcbioRNASeq) <- f



plotGene.DESeqAnalysis <-  # nolint
    function(object) {
        validObject(object)
        # Using DESeqTransform
        dt <- slot(object, "transform")
        if ("rlogIntercept" %in% colnames(dt)) {
            countsAxisLabel <- "rlog counts (log2)"
        } else {
            countsAxisLabel <- "vst counts (log2)"
        }
        do.call(
            what = plotGene,
            args = matchArgsToDoCall(
                args = list(
                    object = dt,
                    genes = genes,
                    countsAxisLabel = countsAxisLabel
                )
            )
        )
    }
f <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f <- f[setdiff(names(f), c("assay", "countsAxisLabel"))]
formals(plotGene.DESeqAnalysis) <- f



#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("bcbioRNASeq"),
    definition = plotGene.bcbioRNASeq
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("DESeqAnalysis"),
    definition = plotGene.DESeqAnalysis
)
