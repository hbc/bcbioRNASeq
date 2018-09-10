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
#'     object = deseq_small@results[[1L]],
#'     counts = deseq_small@transform
#' )
#'
#' # DESeqResults, bcbioRNASeq ====
#' plotDEGHeatmap(
#'     object = deseq_small@results[[1L]],
#'     counts = bcb_small,
#'     normalized = "vst"
#' )
NULL



.plotDEGHeatmap.DESeqResults <-  # nolint
    function(
        object,
        counts,
        scale = "row",
        title = TRUE
        # Defining additional formals below.
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        # Coerce to RSE then SE to preserve rowData.
        if (is(counts, "RangedSummarizedExperiment")) {
            counts <- as(counts, "RangedSummarizedExperiment")
        }
        counts <- as(counts, "SummarizedExperiment")
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
                "counts",
                "alpha",
                "lfcThreshold"
            ),
            call = matchS4Call()
        )
        do.call(what = plotHeatmap, args = args)
    }

# Define the formals.
f1 <- formals(.plotDEGHeatmap.DESeqResults)
f2 <- methodFormals(
    f = "plotHeatmap",
    signature = "SummarizedExperiment"
)
f2 <- f2[setdiff(
    x = names(f2),
    y = c(names(f1), "object")
)]
f <- c(f1, f2)
formals(.plotDEGHeatmap.DESeqResults) <- f



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature(
        object = "DESeqResults",
        counts = "DESeqTransform"
    ),
    definition = .plotDEGHeatmap.DESeqResults
)



# FIXME Improve the formals here.
#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature(
        object = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    definition = function(
        object,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        # Coerce to `DESeqTransform`.
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts(counts, normalized = normalized))
        args <- list(
            object = object,
            counts = DESeqTransform(rse),
            ...
        )
        do.call(what = plotDEGHeatmap, args = args)
    }
)
