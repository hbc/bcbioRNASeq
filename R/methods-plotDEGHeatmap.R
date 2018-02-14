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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#' load(system.file(
#'     file.path("extdata", "dds.rda"),
#'     package = "bcbioRNASeq"))
#' load(system.file(
#'     file.path("extdata", "res.rda"),
#'     package = "bcbioRNASeq"))
#' load(system.file(
#'     file.path("extdata", "rld.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # Use our stashed gene2symbol
#' gene2symbol <- gene2symbol(bcb)
#' annotationCol <- sampleMetadata(bcb) %>%
#'     .[, interestingGroups(bcb), drop = FALSE]
#'
#' # DESeqResults, DESeqTransform
#' plotDEGHeatmap(
#'     res,
#'     counts = rld,
#'     gene2symbol = gene2symbol,
#'     annotationCol = annotationCol)
#'
#' # DESeqResults, DESeqDataSet
#' # Using default ggplot2 colors
#' plotDEGHeatmap(
#'     res,
#'     counts = dds,
#'     gene2symbol = gene2symbol,
#'     color = NULL,
#'     legendColor = NULL)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase camel
.plotDEGHeatmap <- function(
    object,
    counts,
    alpha = 0.01,
    lfc = 0L,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis::viridis,
    legendColor = viridis::viridis,
    title,
    ...) {
    # Passthrough: color, legendColor
    assert_is_data.frame(object)
    assert_are_identical(rownames(object), rownames(counts))
    assert_is_numeric(alpha)
    assert_is_scalar(alpha)
    assert_is_implicit_integer(lfc)
    assert_is_scalar(lfc)
    assert_formal_gene2symbol(object, genes, gene2symbol)
    assert_formal_color_function(color)
    assert_formal_color_function(legendColor)

    results <- object %>%
        camel(strict = FALSE) %>%
        # Keep genes that pass alpha cutoff
        .[!is.na(.[["padj"]]), , drop = FALSE] %>%
        .[.[["padj"]] < alpha, , drop = FALSE] %>%
        # Keep genes that pass log2 fold change cutoff
        .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
        .[.[["log2FoldChange"]] > lfc |
            .[["log2FoldChange"]] < -lfc, , drop = FALSE]

    if (nrow(results) == 0L) {
        warn("No genes passed significance cutoffs")
        return(NULL)
    }

    genes <- rownames(results)
    .checkGenes(genes, gene2symbol)

    counts <- counts[genes, , drop = FALSE]

    .plotHeatmap(
        object = counts,
        gene2symbol = gene2symbol,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...)
}



.plotDEGHeatmap.DESeqResults <- function(  # nolint
    object,
    counts,
    alpha,
    lfc = 0L,
    gene2symbol = TRUE,
    annotationCol = NULL,
    scale = "row",
    color = viridis::viridis,
    legendColor = viridis::viridis,
    title = TRUE,
    ...) {
    results <- as.data.frame(object)
    if (is(counts, "DESeqDataSet") |
        is(counts, "DESeqTransform")) {
        counts <- assay(counts)
    }
    if (missing(alpha)) {
        alpha <- metadata(object)[["alpha"]]
    }
    if (isTRUE(title)) {
        title <- .resContrastName(object)
    }
    .plotDEGHeatmap(
        object = results,
        counts = counts,
        alpha = alpha,
        lfc = lfc,
        gene2symbol = gene2symbol,
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
