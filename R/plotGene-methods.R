#' Plot Gene Expression
#'
#' @name plotGene
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#' @importFrom basejump plotGene
#' @export
#'
#' @inherit basejump::plotGene
#'
#' @inheritParams general
#'
#' @examples
#' g2s <- gene2symbol(bcb_small)
#' geneIDs <- head(g2s[["geneID"]])
#' print(geneIDs)
#' geneNames <- head(g2s[["geneName"]])
#' print(geneNames)
#'
#' # bcbioRNASeq ====
#' object <- bcb_small
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
#' # DESeqAnalysis ====
#' object <- deseq_small
#' plotGene(object, genes = geneIDs, style = "facet")
#' plotGene(object, genes = geneNames, style = "wide")
NULL



.plotGene.bcbioRNASeq <-  # nolint
    function(
        object,
        genes,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle")
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
f1 <- formals(.plotGene.bcbioRNASeq)
f2 <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f2 <- f2[setdiff(names(f2), c(names(f1), "countsAxisLabel"))]
f <- c(f1, f2)
formals(.plotGene.bcbioRNASeq) <- f



.plotGene.DESeqTransform <-  # nolint
    function(object) {
        validObject(object)
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            normalized <- "rlog"
        } else {
            normalized <- "vst"
        }
        do.call(
            what = plotGene,
            args = matchArgsToDoCall(
                args = list(
                    object = as(object, "RangedSummarizedExperiment"),
                    genes = genes,
                    countsAxisLabel = paste(normalized, "counts (log2)")
                )
            )
        )
    }
f <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f <- f[setdiff(names(f), "countsAxisLabel")]
formals(.plotGene.DESeqTransform) <- f



.plotGene.DESeqAnalysis <-  # nolint
    function(object) {
        validObject(object)
        do.call(
            what = plotGene,
            args = matchArgsToDoCall(
                args = list(
                    object = slot(object, "transform"),
                    genes = genes
                )
            )
        )
    }
f <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f <- f[setdiff(names(f), "countsAxisLabel")]
formals(.plotGene.DESeqAnalysis) <- f



#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("bcbioRNASeq"),
    definition = .plotGene.bcbioRNASeq
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("DESeqAnalysis"),
    definition = .plotGene.DESeqAnalysis
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("DESeqTransform"),
    definition = .plotGene.DESeqTransform
)
