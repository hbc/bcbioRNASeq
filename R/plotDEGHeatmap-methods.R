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
    alpha,
    lfcThreshold = 0L,
    scale = c("row", "column", "none"),
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
    if (missing(alpha)) {
        alpha <- metadata(results)[["alpha"]]
    }
    assert_is_a_number(alpha)
    assert_is_a_number(lfcThreshold)
    assert_all_are_non_negative(lfcThreshold)
    scale <- match.arg(scale)

    # Title
    if (isTRUE(title)) {
        title <- contrastName(results)
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    deg <- significants(results, padj = alpha, fc = lfcThreshold)

    # Early return if there are no DEGs.
    if (!length(deg)) {
        warning("No significant DEGs to plot")
        return(invisible())
    }

    # Subset the counts to only contain DEGs.
    counts <- counts[deg, , drop = FALSE]

    # Using SummarizedExperiment method here.
    args <- as.list(matchS4Call())[-1L]
    args[["object"]] <- counts
    args <- args[setdiff(
        x = names(args),
        y = c(
            "results",
            "counts",
            "alpha",
            "lfcThreshold",
            "..."
        )
    )]
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
formals(.plotDEGHeatmap) <- as.pairlist(f)



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
        plotDEGHeatmap(
            results = results,
            counts = rse,
            ...
        )
    }
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
        # FIXME
        print(sys.calls())
        args <- as.list(matchS4Call())[-1L]
        args[["counts"]] <- rse
        do.call(
            what = plotDEGHeatmap,
            args = args
        )
    }
)
