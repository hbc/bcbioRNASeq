# FIXME Use `clusteringMethod = "ward.D2" by default`.
# FIXME Don't message about gene symbol mappings.



# Do not allow post hoc alpha, lfcThreshold cutoffs.



#' Differentially Expressed Gene Heatmap
#' @name plotDEGHeatmap
#' @inherit basejump::plotHeatmap
#' @author Michael Steinbaugh
#'
#' @details
#' This function is a simplified version of [basejump::plotHeatmap()] that is
#' optimized for handling a `DESeqResults` object rather a gene vector. All of
#' the optional parameters for `plotHeatmap` are also available to this
#' function.
#'
#' To adjust the annotation columns, modify the `colData` of the `counts`
#' argument, which must contain a `SummarizedExperiment` (e.g. `DESeqTransform`,
#' `DESeqDataSet`).
#'
#' @inheritParams general
#' @param counts `DESeqTransform`.
#'
#' @examples
#' data(deseq_small)
#'
#' ## DESeqAnalysis ====
#' plotDEGHeatmap(deseq_small)
#'
#' ## DESeqResults ====
#' plotDEGHeatmap(
#'     object = as(deseq_small, "DESeqResults"),
#'     counts = as(deseq_small, "DESeqTransform")
#' )
NULL



# DESeqResults =================================================================
plotDEGHeatmap.DESeqResults <-  # nolint
    function(
        object,
        counts,
        direction = c("both", "up", "down"),
        scale = "row",
        title = TRUE
    ) {
        assert_is_all_of(object, "DESeqResults")
        validObject(object)
        assert_is_all_of(counts, "DESeqTransform")
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        alpha <- metadata(object)[["alpha"]]
        assertIsAlpha(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        direction <- match.arg(direction)
        # Hiding the choices from the user by default, because in most cases row
        # scaling should be used.
        scale <- match.arg(scale, choices = c("row", "column", "none"))

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

        # Coerce DESeqTransform to RSE then SE to preserve rowData.
        if (is(counts, "RangedSummarizedExperiment")) {
            counts <- as(counts, "RangedSummarizedExperiment")
        }
        counts <- as(counts, "SummarizedExperiment")

        # Subset the counts to only contain DEGs.
        counts <- counts[deg, , drop = FALSE]

        # Title
        if (isTRUE(title)) {
            title <- paste0(
                contrastName(object), "\n",
                "(",
                "alpha < ", alpha, "; ",
                "n = ", length(deg),
                ")"
            )
        } else if (!is_a_string(title)) {
            title <- NULL
        }

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
                    "lfcThreshold",
                    "direction"
                )
            )
        )
    }
f1 <- formals(plotDEGHeatmap.DESeqResults)
f2 <- methodFormals(f = "plotHeatmap", signature = "SummarizedExperiment")
f2 <- f2[setdiff(names(f2), c(names(f1), "object", "assay"))]
f <- c(f1, f2)
formals(plotDEGHeatmap.DESeqResults) <- f



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature("DESeqResults"),
    definition = plotDEGHeatmap.DESeqResults
)



# DESeqAnalysis ================================================================
plotDEGHeatmap.DESeqAnalysis <-  # nolint
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
f1 <- formals(plotDEGHeatmap.DESeqAnalysis)
f2 <- formals(plotDEGHeatmap.DESeqResults)
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(plotDEGHeatmap.DESeqAnalysis) <- f



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature("DESeqAnalysis"),
    definition = plotDEGHeatmap.DESeqAnalysis
)
