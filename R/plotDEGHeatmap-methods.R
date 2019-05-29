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
#' @inheritParams general
#' @param ... Passthrough arguments to [plotHeatmap()].
#'
#' @seealso
#' - `help("plotHeatmap", "bcbioBase")`.
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



#' @rdname plotDEGHeatmap
#' @export
setGeneric(
    "plotDEGHeatmap",
    function(results, counts, ...) {
        standardGeneric("plotDEGHeatmap")
    }
)



`plotDEGHeatmap.DESeqResults,SummarizedExperiment` <-  # nolint
    function(
        results,
        counts,
        alpha,
        lfcThreshold = 0L,
        title = TRUE,
        ...
    ) {
        validObject(results)
        validObject(counts)
        # Coerce to RSE then SE to preserve rowData
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

        # Title
        if (isTRUE(title)) {
            title <- contrastName(results)
        } else if (!is_a_string(title)) {
            title <- NULL
        }

        deg <- significants(results, padj = alpha, fc = lfcThreshold)

        # Early return if there are no DEGs
        if (!length(deg)) {
            warning("No significant DEGs to plot")
            return(invisible())
        }

        # Subset the counts to only contain DEGs
        counts <- counts[deg, , drop = FALSE]

        # SummarizedExperiment method
        plotHeatmap(
            object = counts,
            title = title,
            ...
        )
    }



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature(
        results = "DESeqResults",
        counts = "SummarizedExperiment"
    ),
    definition = `plotDEGHeatmap.DESeqResults,SummarizedExperiment`
)



`plotDEGHeatmap.DESeqResults,bcbioRNASeq` <-  # nolint
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



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature(
        results = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    definition = `plotDEGHeatmap.DESeqResults,bcbioRNASeq`
)



`plotDEGHeatmap.DESeqResults,DESeqDataSet` <-  # nolint
    function(
        results,
        counts,
        ...
    ) {
        validObject(counts)
        message("Using normalized counts")
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts(counts, normalized = TRUE))
        plotDEGHeatmap(
            results = results,
            counts = rse,
            ...
        )
    }



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature(
        results = "DESeqResults",
        counts = "DESeqDataSet"
    ),
    definition = `plotDEGHeatmap.DESeqResults,DESeqDataSet`
)



`plotDEGHeatmap.DESeqResults,DESeqTransform` <-  # nolint
    `plotDEGHeatmap.DESeqResults,SummarizedExperiment`



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    f = "plotDEGHeatmap",
    signature = signature(
        results = "DESeqResults",
        counts = "DESeqTransform"
    ),
    definition = `plotDEGHeatmap.DESeqResults,DESeqTransform`
)
