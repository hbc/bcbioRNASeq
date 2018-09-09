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
#'
#' @seealso
#' - `help("plotHeatmap", "basejump")`.
#' - `findMethod("plotHeatmap", "SummarizedExperiment")`.
#'
#' @examples
#' # DESeqResults, bcbioRNASeq ====
#' plotDEGHeatmap(
#'     results = res_small,
#'     counts = bcb_small,
#'     normalized = "vst"
#' )
#'
#' # DESeqResults, DESeqTransform ====
#' vst_small <- DESeq2::varianceStabilizingTransformation(dds_small)
#' plotDEGHeatmap(
#'     results = res_small,
#'     counts = vst_small
#' )
#'
#' # DESeqResults, DESeqDataSet ====
#' # This always uses normalized counts
#' # Using default ggplot2 colors
#' plotDEGHeatmap(
#'     results = res_small,
#'     counts = dds_small,
#'     color = NULL,
#'     legendColor = NULL
#' )
NULL



.plotDEGHeatmap <-  # nolint
function(
    results,
    counts,
    alpha = NULL,
    lfcThreshold = 0L,
    scale = "row",
    title = TRUE
    # Defining additional formals below.
) {
    validObject(results)
    validObject(counts)
    # Coerce to RSE then SE to preserve rowData.
    if (is(counts, "RangedSummarizedExperiment")) {
        counts <- as(counts, "RangedSummarizedExperiment")
    }
    counts <- as(counts, "SummarizedExperiment")
    assert_are_identical(rownames(results), rownames(counts))
    if (is.null(alpha)) {
        alpha <- metadata(results)[["alpha"]]
    }
    assert_is_a_number(alpha)
    assert_is_a_number(lfcThreshold)
    assert_all_are_non_negative(lfcThreshold)
    # Hiding the choices from the user by default, because in most cases row
    # scaling should be used.
    scale <- match.arg(
        arg = scale,
        choices = c("row", "column", "none")
    )

    # Title
    if (isTRUE(title)) {
        title <- contrastName(results)
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    deg <- significants(results, padj = alpha, fc = lfcThreshold)

    # Early return if there are no DEGs.
    if (!length(deg)) {
        warning("No significant DEGs to plot", call. = FALSE)
        return(invisible())
    }

    # Subset the counts to only contain DEGs.
    counts <- counts[deg, , drop = FALSE]

    # Using `do.call()` return with SummarizedExperiment method here.
    args <- setArgsToDoCall(
        args = list(
            object = counts,
            interestingGroups = interestingGroups,
            scale = scale,
            title = title,
            ...
        ),
        removeArgs = c(
            "results",
            "counts",
            "alpha",
            "lfcThreshold"
        ),
        call = matchS4Call(),
        fun = sys.function()
    )
    do.call(what = plotHeatmap, args = args)
}

# Define the formals.
f1 <- formals(.plotDEGHeatmap)
f2 <- methodFormals(
    f = "plotHeatmap",
    signature = "SummarizedExperiment"
)
f2 <- f2[setdiff(
    x = names(f2),
    y = c(names(f1), "object")
)]
f <- c(f1, f2)
formals(.plotDEGHeatmap) <- f



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        results = "DESeqResults",
        counts = "SummarizedExperiment"
    ),
    .plotDEGHeatmap
)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        results = "DESeqResults",
        counts = "DESeqTransform"
    ),
    getMethod(
        "plotDEGHeatmap",
        signature(
            results = "DESeqResults",
            counts = "SummarizedExperiment"
        )
    )
)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        results = "DESeqResults",
        counts = "DESeqDataSet"
    ),
    function(
        results,
        counts,
        ...
    ) {
        validObject(counts)
        message("Using normalized counts")
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts(counts, normalized = TRUE))
        args <- list(
            results = results,
            counts = rse,
            ...
        )
        do.call(what = plotDEGHeatmap, args = args)
    }
)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        results = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    function(
        results,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(counts)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts(counts, normalized = normalized))
        args <- list(
            results = results,
            counts = rse,
            ...
        )
        do.call(what = plotDEGHeatmap, args = args)
    }
)
