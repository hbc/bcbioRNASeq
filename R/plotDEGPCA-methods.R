# TODO Show lfcThreshold info on the plot.



#' @name plotDEGPCA
#' @inherit basejump::plotDEGPCA
#' @inherit plotPCA
#' @include plotPCA-methods.R
#' @author Michael Steinbaugh
#'
#' @inheritParams params
#' @inheritParams basejump::params
#' @param counts `DESeqTransform`.
#'
#' @examples
#' data(deseq)
#'
#' ## DESeqAnalysis ====
#' plotDEGPCA(deseq)
#'
#' ## DESeqResults ====
#' plotDEGPCA(
#'     object = as(deseq, "DESeqResults"),
#'     counts = as(deseq, "DESeqTransform")
#' )
NULL



#' @importFrom basejump plotDEGPCA
#' @aliases NULL
#' @export
basejump::plotDEGPCA



# DESeqResults =================================================================
# Do not allow post hoc alpha, lfcThreshold cutoffs.
plotDEGPCA.DESeqResults <-  # nolint
    function(
        object,
        counts,
        direction = c("both", "up", "down")
    ) {
        assert_is_all_of(object, "DESeqResults")
        validObject(object)
        assert_is_all_of(counts, "DESeqTransform")
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        interestingGroups <- matchInterestingGroups(
            object = counts,
            interestingGroups = interestingGroups
        )
        interestingGroups(counts) <- interestingGroups
        alpha <- metadata(object)[["alpha"]]
        assertIsAlpha(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        assertIsColorScaleDiscreteOrNULL(color)
        direction <- match.arg(direction)
        assert_is_a_bool(label)
        return <- match.arg(return)

        # Get the character vector of DEGs.
        deg <- .deg(
            object = object,
            alpha = alpha,
            lfcThreshold = lfcThreshold,
            direction = direction
        )
        if (!has_length(deg)) {
            warning("No significant DEGs to plot.", call. = FALSE)
            return(invisible())
        }

        # Using our internal `SummarizedExperiment` method.
        do.call(
            what = plotPCA.SummarizedExperiment,
            args = list(
                object = counts[deg, , drop = FALSE],
                interestingGroups = interestingGroups,
                ntop = Inf,
                label = label,
                title = contrastName(object),
                subtitle = paste(
                    length(deg), "genes;",
                    "alpha <", alpha
                ),
                return = return
            )
        )
    }
f1 <- formals(plotDEGPCA.DESeqResults)
# Note that we're not exporting the plotPCA SE method.
f2 <- formals(plotPCA.SummarizedExperiment)
f2 <- f2[c("interestingGroups", "color", "label", "return")]
f <- c(f1, f2)
formals(plotDEGPCA.DESeqResults) <- f



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature("DESeqResults"),
    definition = plotDEGPCA.DESeqResults
)



# DESeqAnalysis ================================================================
plotDEGPCA.DESeqAnalysis <-  # nolint
    function(
        object,
        counts = NULL,
        results
    ) {
        results <- .matchResults(
            object = object,
            results = results
        )
        counts <- slot(object, "transform")
        do.call(
            what = plotDEGPCA,
            args = matchArgsToDoCall(
                args = list(
                    object = results,
                    counts = counts
                ),
                removeFormals = "results"
            )
        )
    }
f1 <- formals(plotDEGPCA.DESeqAnalysis)
f2 <- formals(plotDEGPCA.DESeqResults)
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(plotDEGPCA.DESeqAnalysis) <- f



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature("DESeqAnalysis"),
    definition = plotDEGPCA.DESeqAnalysis
)
