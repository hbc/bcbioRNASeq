#' @name plotCountsPerGene
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @importMethodsFrom minimalism plotCountsPerGene
#' @inherit bioverbs::plotCountsPerGene
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @section Trimmed Mean of M-Values:
#'
#' We recommend visualizing counts normalized with the **T**rimmed **M**ean of
#' **M**-Values (TMM) method here. TMM normalization equates the overall
#' expression levels of genes between samples under the assumption that the
#' majority of them are not differentially expressed. Therefore, by normalizing
#' for total RNA expression by sample, we expect the spread of the
#' TMM-normalized counts per gene to be similar for every sample.
#'
#' @references TMM: Robinson, et al., 2010.
#'
#' @examples
#' data(bcb)
#' plotCountsPerGene(bcb)
NULL



#' @importFrom bioverbs plotCountsPerGene
#' @aliases NULL
#' @export
bioverbs::plotCountsPerGene



plotCountsPerGene.bcbioRNASeq <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)

        # Coerce to RSE.
        rse <- as(object, "RangedSummarizedExperiment")
        counts <- counts(object, normalized = normalized)
        assays(rse) <- list(counts)
        assayNames(rse) <- normalized

        # Set the counts axis label.
        countsAxisLabel <- paste(normalized, "counts")
        trans <- .normalizedTrans(normalized)
        if (trans != "identity") {
            countsAxisLabel <- paste(trans, countsAxisLabel)
        }

        do.call(
            what = plotCountsPerGene,
            args = matchArgsToDoCall(
                args = list(
                    object = rse,
                    trans = trans,
                    countsAxisLabel = countsAxisLabel
                ),
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(plotCountsPerGene.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCountsPerGene",
    signature = "SummarizedExperiment",
    package = "minimalism"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
# Ensure TMM is set first.
f[["normalized"]] <- unique(c("tmm", normalizedCounts))
formals(plotCountsPerGene.bcbioRNASeq) <- f



#' @rdname plotCountsPerGene
#' @export
setMethod(
    f = "plotCountsPerGene",
    signature = signature("bcbioRNASeq"),
    definition = plotCountsPerGene.bcbioRNASeq
)
