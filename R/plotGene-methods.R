# FIXME Need to improve the formals.
# FIXME Add `DESeqAnalysis` support.
# FIXME Add DESeqAnalysis method support.



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
#' @examples
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
#'     return = "facet"
#' )
#' plotGene(
#'     object = object,
#'     genes = geneNames,
#'     normalized = "vst",
#'     return = "wide"
#' )
#'
#' # DESeqAnalysis ====
#' object <- deseq_small
#' plotGene(object, genes = geneIDs, return = "facet")
#' plotGene(object, genes = geneNames, return = "wide")
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
        args <- setArgsToDoCall(
            args = list(
                object = rse,
                genes = genes,
                countsAxisLabel = paste(normalized, "counts (log2)")
            ),
            removeArgs = "normalized",
            call = matchCall()
        )
        do.call(what = plotGene, args = args)
    }

# Assign the formals.
f1 <- formals(.plotGene.bcbioRNASeq)
f2 <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f2 <- f2[setdiff(names(f2), c(names(f1), "countsAxisLabel"))]
f <- c(f1, f2)
formals(.plotGene.bcbioRNASeq) <- f



.plotGene.DESeqDataSet <-  # nolint
    function() {
        validObject(object)
        counts <- counts(object, normalized = TRUE)
        # Ensure counts are log2 scale.
        counts <- log2(counts + 1L)
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts)
        args <- setArgsToDoCall(
            args = list(
                object = rse,
                genes = genes,
                countsAxisLabel = "normalized counts (log2)"
            ),
            removeArgs = "normalized",
            call = matchCall()
        )
        do.call(what = plotGene, args = args)
    }

# Assign the formals.
f <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f <- f[setdiff(names(f), "countsAxisLabel")]
formals(.plotGene.DESeqDataSet) <- f



.plotGene.DESeqTransform <-  # nolint
    function() {
        validObject(object)
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            normalized <- "rlog"
        } else {
            normalized <- "vst"
        }
        args <- setArgsToDoCall(
            args = list(
                object = as(object, "RangedSummarizedExperiment"),
                genes = genes,
                countsAxisLabel = paste(normalized, "counts (log2)")
            ),
            call = matchCall()
        )
        do.call(what = plotGene, args = args)
    }

# Assign the formals.
f <- methodFormals(f = "plotGene", signature = "SummarizedExperiment")
f <- f[setdiff(names(f), "countsAxisLabel")]
formals(.plotGene.DESeqTransform) <- f



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
    signature = signature("DESeqTransform"),
    definition = .plotGene.DESeqTransform
)



#' @rdname plotGene
#' @export
setMethod(
    f = "plotGene",
    signature = signature("DESeqDataSet"),
    definition = .plotGene.DESeqDataSet
)
