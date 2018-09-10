# FIXME Define and recommend DESeqAnalysis method
# FIXME Include the alpha on the plot



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
#' # DESeqResults, DESeqTransform ====
#' plotDEGHeatmap(
#'     object = deseq_small@results,
#'     counts = deseq_small@transform
#' )
NULL



.plotDEGHeatmap <-  # nolint
function(
    object,
    counts,
    alpha = NULL,
    lfcThreshold = 0L,
    scale = "row",
    title = TRUE
    # Defining additional formals below.
) {
    validObject(object)
    validObject(counts)
    assert_are_identical(
        x = rownames(object),
        y = rownames(counts)
    )
    # Coerce to RSE then SE to preserve rowData.
    if (is(counts, "RangedSummarizedExperiment")) {
        counts <- as(counts, "RangedSummarizedExperiment")
    }
    counts <- as(counts, "SummarizedExperiment")
    if (is.null(alpha)) {
        alpha <- metadata(object)[["alpha"]]
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
        title <- contrastName(object)
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    deg <- significants(object, padj = alpha, fc = lfcThreshold)

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
            "object",
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
        object = "DESeqResults",
        counts = "SummarizedExperiment"
    ),
    .plotDEGHeatmap
)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        object = "DESeqResults",
        counts = "DESeqTransform"
    ),
    getMethod(
        "plotDEGHeatmap",
        signature(
            object = "DESeqResults",
            counts = "SummarizedExperiment"
        )
    )
)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        object = "DESeqResults",
        counts = "DESeqDataSet"
    ),
    function(
        object,
        counts,
        ...
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(
            x = rownames(object),
            y = rownames(counts)
        )
        message("Using normalized counts")
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts(counts, normalized = TRUE))
        args <- list(
            object = object,
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
        object = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    function(
        object,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(
            x = rownames(object),
            y = rownames(counts)
        )
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts(counts, normalized = normalized))
        args <- list(
            object = object,
            counts = rse,
            ...
        )
        do.call(what = plotDEGHeatmap, args = args)
    }
)
