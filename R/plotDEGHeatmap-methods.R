#' Differentially Expressed Gene Heatmap
#'
#' This function is a simplified version of [plotHeatmap()] that is
#' optimized for handling a `DESeqResults` object rather a gene vector. All of
#' the optional parameters for [plotHeatmap()] are also available to this
#' function.
#'
#' To adjust the annotation columns, modify the `colData` of the `counts`
#' argument, which must contain a `SummarizedExperiment` (e.g. `DESeqTransform`,
#' `DESeqDataSet`).
#'
#' @name plotDEGHeatmap
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inherit basejump::plotHeatmap
#'
#' @inheritParams general
#' @param counts `DESeqTransform`.
#'
#' @seealso
#' - `help("plotHeatmap", "basejump")`.
#' - `findMethod("plotHeatmap", "SummarizedExperiment")`.
#'
#' @examples
#' # DESeqAnalysis ====
#' plotDEGHeatmap(deseq_small)
#'
#' # DESeqResults ====
#' plotDEGHeatmap(
#'     object = as(deseq_small, "DESeqResults"),
#'     counts = as(deseq_small, "DESeqTransform")
#' )
NULL



.plotDEGHeatmap.DESeqResults <-  # nolint
    function(
        object,
        counts,
        scale = "row",
        title = TRUE
    ) {
        assert_is_all_of(object, "DESeqResults")
        validObject(object)
        assert_is_all_of(counts, "DESeqTransform")
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        alpha <- metadata(object)[["alpha"]]
        assert_is_a_number(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        # Hiding the choices from the user by default, because in most cases row
        # scaling should be used.
        scale <- match.arg(scale, choices = c("row", "column", "none"))

        # Title
        if (isTRUE(title)) {
            title <- paste0(contrastName(object), " (alpha < ", alpha, ")")
        } else if (!is_a_string(title)) {
            title <- NULL
        }

        deg <- significants(object, padj = alpha, fc = lfcThreshold)

        # Early return if there are no DEGs.
        if (!length(deg) > 0L) {
            warning("No significant DEGs to plot", call. = FALSE)
            return(invisible())
        }

        # Coerce DESeqTransform to RSE then SE to preserve rowData.
        if (is(counts, "RangedSummarizedExperiment")) {
            counts <- as(counts, "RangedSummarizedExperiment")
        }
        counts <- as(counts, "SummarizedExperiment")

        # Subset the counts to only contain DEGs.
        counts <- counts[deg, , drop = FALSE]

        # Using `do.call()` return with SummarizedExperiment method here.
        do.call(
            what = plotHeatmap,
            args = matchArgsToDoCall(
                args = list(
                    object = counts,
                    title = title
                ),
                removeFormals = c(
                    "counts",
                    "alpha",
                    "lfcThreshold"
                )
            )
        )
    }
f1 <- formals(.plotDEGHeatmap.DESeqResults)
f2 <- methodFormals(f = "plotHeatmap", signature = "SummarizedExperiment")
f2 <- f2[setdiff(names(f2), c(names(f1), "object"))]
f <- c(f1, f2)
formals(.plotDEGHeatmap.DESeqResults) <- f



.plotDEGHeatmap.DESeqAnalysis <-  # nolint
    function(object, results) {
        do.call(
            what = plotDEGHeatmap,
            args = matchArgsToDoCall(
                args = list(
                    # DESeqResults
                    object = .matchResults(
                        object = object,
                        results = results
                    ),
                    # DESeqTransform
                    counts = as(object, "DESeqTransform")
                ),
                removeFormals = "results"
            )
        )
    }
f1 <- formals(.plotDEGHeatmap.DESeqAnalysis)
f2 <- formals(.plotDEGHeatmap.DESeqResults)
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(.plotDEGHeatmap.DESeqAnalysis) <- f



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature("DESeqAnalysis"),
    definition = .plotDEGHeatmap.DESeqAnalysis
)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature("DESeqResults"),
    definition = .plotDEGHeatmap.DESeqResults
)
