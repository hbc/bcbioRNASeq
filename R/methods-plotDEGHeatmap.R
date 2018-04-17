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
#' @inherit bcbioBase::plotHeatmap
#'
#' @inheritParams general
#' @param ... Passthrough arguments to [plotHeatmap()].
#'
#' @seealso
#' - `help("plotHeatmap", "bcbioBase")`.
#' - `findMethod("plotHeatmap", "SummarizedExperiment")`.
#'
#' @examples
#' # DESeqResults, SummarizedExperiment ====
#' plotDEGHeatmap(
#'     results = res_small,
#'     counts = rld_small
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
#'
#' # DESeqResults, bcbioRNASeq ====
#' plotDEGHeatmap(
#'     results = res_small,
#'     counts = bcb_small
#' )
NULL



# Methods ======================================================================
#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        results = "DESeqResults",
        counts = "SummarizedExperiment"
    ),
    function(
        results,
        counts,
        lfc = 0L,
        title = TRUE,
        ...
    ) {
        # Legacy arguments =====================================================
        call <- match.call(expand.dots = TRUE)
        # gene2symbol
        if ("genesymbol" %in% names(call)) {
            stop(paste(
                "`gene2symbol` argument has been deprecated in favor of",
                "stashing gene-to-symbol mappings in `rowRanges`"
            ))
        }

        # Assert checks ========================================================
        validObject(results)
        assert_are_identical(rownames(results), rownames(counts))
        assert_is_a_number(lfc)
        assert_all_are_non_negative(lfc)

        # Title
        if (isTRUE(title)) {
            title <- contrastName(results)
        } else if (!is_a_string(title)) {
            title <- NULL
        }

        # Alpha level
        alpha <- metadata(results)[["alpha"]]
        assert_is_a_number(alpha)

        # data.frame containing only the significant DEGs
        data <- results %>%
            as.data.frame() %>%
            camel() %>%
            # Keep genes that pass alpha cutoff
            .[!is.na(.[["padj"]]), , drop = FALSE] %>%
            .[.[["padj"]] < alpha, , drop = FALSE] %>%
            # Keep genes that pass log2 fold change cutoff
            .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
            .[.[["log2FoldChange"]] > lfc |
                  .[["log2FoldChange"]] < -lfc, , drop = FALSE]
        assert_has_rows(data)
        deg <- rownames(data)

        # Subset the counts to only contain DEGs
        counts <- counts[deg, , drop = FALSE]

        # SummarizedExperiment method
        plotHeatmap(
            object = counts,
            title = title,
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
        assay(rse) <- counts(counts, normalized = TRUE)
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
        counts = "bcbioRNASeq"
    ),
    function(
        results,
        counts,
        normalized = c("rlog", "vst", "tmm", "tpm"),
        ...
    ) {
        validObject(counts)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(counts, "RangedSummarizedExperiment")
        assay(rse) <- counts(counts, normalized = normalized)
        plotDEGHeatmap(
            results = results,
            counts = rse,
            ...
        )
    }
)
