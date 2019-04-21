#' @name plotCountsPerFeature
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @importMethodsFrom acidplots plotCountsPerFeature
#'
#' @inherit bioverbs::plotCountsPerFeature
#' @inheritParams acidplots::params
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
#' plotCountsPerFeature(bcb)
NULL



#' @rdname plotCountsPerFeature
#' @name plotCountsPerFeature
#' @importFrom bioverbs plotCountsPerFeature
#' @export
NULL



plotCountsPerFeature.bcbioRNASeq <-  # nolint
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
            what = plotCountsPerFeature,
            args = matchArgsToDoCall(
                args = args,
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(plotCountsPerFeature.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCountsPerFeature",
    signature = "SummarizedExperiment",
    package = "acidplots"
)
f2 <- f2[setdiff(
    x = names(f2),
    y = c(names(f1), "assay", "countsAxisLabel", "trans")
)]
f <- c(f1, f2)
# Ensure TPM is set first.
f[["normalized"]] <- unique(c("tmm", normalizedCounts))
f[["trans"]] <- trans
formals(plotCountsPerFeature.bcbioRNASeq) <- f



#' @rdname plotCountsPerFeature
#' @export
setMethod(
    f = "plotCountsPerFeature",
    signature = signature("bcbioRNASeq"),
    definition = plotCountsPerFeature.bcbioRNASeq
)
