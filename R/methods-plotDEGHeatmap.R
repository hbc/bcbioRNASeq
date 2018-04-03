#' Differentially Expressed Gene Heatmap
#'
#' This function is a simplified version of [plotHeatmap()] that is
#' optimized for handling a [DESeqResults] object rather a gene vector. All of
#' the optional parameters for [plotHeatmap()] are also available to this
#' function.
#'
#' @rdname plotDEGHeatmap
#' @name plotDEGHeatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#'
#' @inherit plotHeatmap
#'
#' @param counts Secondary object containing a normalized counts matrix.
#' @param alpha *Optional* Alpha level cutoff. If missing, the function will
#'   use the alpha level defined in the object.
#' @param lfc log2 fold change ratio cutoff.
#' @param ... Passthrough arguments to [plotHeatmap()].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/res.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))
#'
#' # Use our stashed gene2symbol
#' gene2symbol <- gene2symbol(bcb)
#' annotationCol <- sampleMetadata(bcb) %>%
#'     .[, interestingGroups(bcb), drop = FALSE]
#'
#' # DESeqResults, DESeqTransform ====
#' plotDEGHeatmap(
#'     object = res,
#'     counts = rld,
#'     gene2symbol = gene2symbol,
#'     annotationCol = annotationCol)
#'
#' # DESeqResults, DESeqDataSet ====
#' # Using default ggplot2 colors
#' plotDEGHeatmap(
#'     object = res,
#'     counts = dds,
#'     gene2symbol = gene2symbol,
#'     color = NULL,
#'     legendColor = NULL)
NULL



# Constructors =================================================================
.plotDEGHeatmap.DESeqResults <- function(  # nolint
    object,
    counts,
    alpha,
    lfc = 0L,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = c("row", "column", "none"),
    color = viridis,
    legendColor = viridis,
    title = TRUE,
    ...) {
    # Passthrough: annotationCol, color, legendColor
    stopifnot(is(object, "DESeqResults"))
    if (is(counts, "DESeqDataSet") || is(counts, "DESeqTransform")) {
        counts <- assay(counts)
    }
    assert_is_matrix(counts)
    assert_are_identical(rownames(object), rownames(counts))
    if (missing(alpha)) {
        alpha <- metadata(object)[["alpha"]]
    }
    assert_is_a_number(alpha)
    assert_all_are_positive(alpha)
    assertIsAnImplicitInteger(lfc)
    assert_all_are_non_negative(lfc)
    assertFormalGene2symbol(
        x = object,
        genes = rownames(counts),
        gene2symbol = gene2symbol
    )
    scale <- match.arg(scale)
    if (isTRUE(title)) {
        title <- .contrastName.DESeqResults(object)
    } else {
        title <- NULL
    }

    # Subset the counts matrix to match only the significant DEGs
    results <- object %>%
        as.data.frame() %>%
        camel() %>%
        # Keep genes that pass alpha cutoff
        .[!is.na(.[["padj"]]), , drop = FALSE] %>%
        .[.[["padj"]] < alpha, , drop = FALSE] %>%
        # Keep genes that pass log2 fold change cutoff
        .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
        .[.[["log2FoldChange"]] > lfc |
              .[["log2FoldChange"]] < -lfc, , drop = FALSE]
    assert_has_rows(results)
    counts <- counts[rownames(results), , drop = FALSE]

    # Rename rows using gene symbols, if desired
    if (is.data.frame(gene2symbol)) {
        rownames <- gene2symbol[rownames(counts), "symbol", drop = TRUE]
        rownames(counts) <- rownames
    }

    plotHeatmap(
        object = counts,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...)
}



# Methods ======================================================================
#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        object = "DESeqResults",
        counts = "DESeqDataSet"),
    .plotDEGHeatmap.DESeqResults)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        object = "DESeqResults",
        counts = "DESeqTransform"),
    .plotDEGHeatmap.DESeqResults)



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        object = "DESeqResults",
        counts = "matrix"),
    .plotDEGHeatmap.DESeqResults)
