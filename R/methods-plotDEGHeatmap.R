#' Differentially Expressed Gene Heatmap
#'
#' This function is a simplified version of [plotHeatmap()] that is
#' optimized for handling a `DESeqResults` object rather a gene vector. All of
#' the optional parameters for [plotHeatmap()] are also available to this
#' function.
#'
#' @name plotDEGHeatmap
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inherit plotHeatmap
#'
#' @inheritParams general
#' @param ... Passthrough arguments to [plotHeatmap()].
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/res_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/rld_small.rda", package = "bcbioRNASeq"))
#'
#' # Use our stashed gene2symbol
#' gene2symbol <- gene2symbol(bcb_small)
#' annotationCol <- colData(bcb_small) %>%
#'     .[, interestingGroups(bcb_small), drop = FALSE]
#'
#' # DESeqResults, DESeqTransform ====
#' plotDEGHeatmap(
#'     object = res_small,
#'     counts = rld_small,
#'     gene2symbol = gene2symbol,
#'     annotationCol = annotationCol
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



# Constructors =================================================================
#' @importFrom basejump convertGenesToSymbols plotHeatmap
.plotDEGHeatmap <- function(  # nolint
    object,
    counts,
    lfc = 0L,
    gene2symbol = NULL,
    title = TRUE,
    ...
) {
    assert_is_any_of(
        x = counts,
        classes = c("DESeqDataSet", "DESeqTransform", "matrix")
    )
    if (is(counts, "DESeqDataSet") || is(counts, "DESeqTransform")) {
        counts <- assay(counts)
    }
    assert_is_matrix(counts)
    assert_are_identical(rownames(object), rownames(counts))
    assert_is_a_number(lfc)
    assert_all_are_non_negative(lfc)
    assertFormalGene2symbol(object, rownames(counts), gene2symbol)

    # Title
    if (isTRUE(title)) {
        title <- contrastName(object)
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    # Alpha level
    alpha <- metadata(object)[["alpha"]]
    assert_is_a_number(alpha)

    data <- object %>%
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
    counts <- counts[rownames(data), , drop = FALSE]

    # Remap gene identifier rows to symbols
    if (is.data.frame(gene2symbol)) {
        rownames(counts) <- convertGenesToSymbols(
            rownames(counts),
            gene2symbol = gene2symbol
        )
    }

    plotHeatmap(counts, title = title, ...)
}



# Methods ======================================================================
#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature("DESeqResults"),
    .plotDEGHeatmap
)
