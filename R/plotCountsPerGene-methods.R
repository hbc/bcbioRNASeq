#' @name plotCountsPerGene
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @importMethodsFrom minimalism plotCountsPerGene
#'
#' @inherit bioverbs::plotCountsPerGene
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @param trans `character(1)`.
#'   Logarithmic transformation to apply. Note that `vst` and `rlog` counts are
#'   already log2.
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



#' @rdname plotCountsPerGene
#' @name plotCountsPerGene
#' @importFrom bioverbs plotCountsPerGene
#' @export
NULL



plotCountsPerGene.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized,
        trans
    ) {
        normalized <- match.arg(normalized)
        trans <- match.arg(trans)
        args <- .dynamicTrans(
            object = object,
            normalized = normalized,
            trans = trans
        )
        do.call(
            what = plotCountsPerGene,
            args = matchArgsToDoCall(
                args = args,
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
f2 <- f2[setdiff(
    x = names(f2),
    y = c(names(f1), "assay", "countsAxisLabel", "trans")
)]
f <- c(f1, f2)
# Ensure TPM is set first.
f[["normalized"]] <- unique(c("tmm", normalizedCounts))
f[["trans"]] <- trans
formals(plotCountsPerGene.bcbioRNASeq) <- f



#' @rdname plotCountsPerGene
#' @export
setMethod(
    f = "plotCountsPerGene",
    signature = signature("bcbioRNASeq"),
    definition = plotCountsPerGene.bcbioRNASeq
)
