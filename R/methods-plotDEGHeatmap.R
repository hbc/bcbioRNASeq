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
#' @param ... Passthrough arguments to [plotHeatmap()] `SummarizedExperiment`
#'   method.
#'
#' @seealso
#' - `help("plotHeatmap", "bcbioBase")`.
#' - `findMethod("plotHeatmap", "SummarizedExperiment")`.
#'
#' @examples
#' # Use our stashed gene2symbol
#' gene2symbol <- gene2symbol(bcb_small)
#'
#' # DESeqResults, DESeqTransform ====
#' plotDEGHeatmap(
#'     object = res_small,
#'     counts = rld_small,
#'     gene2symbol = gene2symbol
#' )
#'
#' # DESeqResults, DESeqDataSet ====
#' # Using default ggplot2 colors
#' plotDEGHeatmap(
#'     object = res_small,
#'     counts = dds_small,
#'     gene2symbol = gene2symbol,
#'     color = NULL,
#'     legendColor = NULL
#' )
NULL



# Methods ======================================================================
#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature("DESeqResults"),
    function(
        object,
        counts,
        lfc = 0L,
        gene2symbol = NULL,
        title = TRUE,
        ...
    ) {
        validObject(object)
        assert_is_all_of(counts, "SummarizedExperiment")
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        assert_is_a_number(lfc)
        assert_all_are_non_negative(lfc)
        assertFormalGene2symbol(
            x = object,
            genes = rownames(counts),
            gene2symbol = gene2symbol
        )

        # Title
        if (isTRUE(title)) {
            title <- contrastName(object)
        } else if (!is_a_string(title)) {
            title <- NULL
        }

        # Alpha level
        alpha <- metadata(object)[["alpha"]]
        assert_is_a_number(alpha)

        deg <- object %>%
            as.data.frame() %>%
            camel() %>%
            # Keep genes that pass alpha cutoff
            .[!is.na(.[["padj"]]), , drop = FALSE] %>%
            .[.[["padj"]] < alpha, , drop = FALSE] %>%
            # Keep genes that pass log2 fold change cutoff
            .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
            .[.[["log2FoldChange"]] > lfc |
                  .[["log2FoldChange"]] < -lfc, , drop = FALSE]
        assert_has_rows(deg)

        # Subset the counts matrix to only contain DEGs
        counts <- counts[rownames(deg), ]
        validObject(counts)

        # Remap gene identifier rows to symbols
        if (is.data.frame(gene2symbol)) {
            rownames <- convertGenesToSymbols(
                object = rownames(counts),
                gene2symbol = gene2symbol
            )
            rownames(counts) <- rownames
        }

        # SummarizedExperiment method
        plotHeatmap(
            counts,
            title = title,
            ...
        )
    }
)
