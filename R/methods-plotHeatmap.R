#' Plot Heatmap
#'
#' @details
#' When column clustering is disabled, the columns are sorted by the interesting
#' groups (`interestingGroups`) specified in the [bcbioRNASeq] and then the
#' sample names.
#'
#' @rdname plotHeatmap
#' @name plotHeatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotHeatmap
#'
#' @inheritParams general
#' @inheritParams counts
#' @inheritParams gene2symbol
#' @inheritParams plotGene
#'
#' @param samples *Optional.* Samples (colnames) to plot.
#' @param annotationCol *Optional.* [data.frame] that defines annotation
#'   mappings for the columns.
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the row direction or the column direction, or none. Corresponding
#'   values are "row", "column" and "none".
#' @param legendColor Colors to use for legend labels. Defaults to the
#'   [viridis()].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))
#'
#' gene2symbol <- gene2symbol(bcb)
#'
#' # bcbioRNASeq ====
#' plotHeatmap(bcb, genes = head(rownames(bcb), 20L))
#'
#' # Full transcriptome heatmap with default pheatmap colors
#' plotHeatmap(bcb, color = inferno, legendColor = inferno)
#'
#' # DESeqDataSet ====
#' dds <- bcbio(bcb, "DESeqDataSet")
#' plotHeatmap(
#'     dds,
#'     genes = head(rownames(dds), 20L),
#'     gene2symbol = gene2symbol)
#'
#' # DESeqTransform ====
#' rld <- assays(bcb)[["rlog"]]
#' plotHeatmap(
#'     rld,
#'     genes = head(rownames(dds), 20L),
#'     gene2symbol = gene2symbol)
NULL



# Constructors =================================================================
.plotHeatmap <- function(
    object,
    samples = NULL,
    genes = NULL,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis,
    legendColor = viridis,
    title = NULL,
    ...) {
    assert_is_matrix(object)
    assert_is_character_or_null(samples)
    assert_is_character_or_null(genes)
    assert_formal_gene2symbol(object, genes, gene2symbol)
    assert_formal_annotation_col(object, annotationCol)
    assert_formal_color_function(color)
    assert_formal_color_function(legendColor)
    assert_is_a_string_or_null(title)

    # Resize the counts matrix
    if (is.vector(samples)) {
        object <- object[, samples, drop = FALSE]
    }
    if (is.vector(genes)) {
        object <- object[genes, , drop = FALSE]
    }

    # Set the rownames to gene symbols
    if (is.data.frame(gene2symbol)) {
        rownames(object) <- rownames(object) %>%
            gene2symbol[., "symbol", drop = TRUE] %>%
            make.names(unique = TRUE)
    }

    plotHeatmap(
        object = object,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...)
}



.plotHeatmap.bcbioRNASeq <- function(  # nolint
    object,
    normalized = "rlog",
    samples = NULL,
    genes = NULL,
    scale = "row",
    color = viridis,
    legendColor = viridis,
    title = NULL,
    ...) {
    counts <- counts(object, normalized = normalized)
    annotationCol <- colData(object) %>%
        .[colnames(counts), interestingGroups(object), drop = FALSE] %>%
        as.data.frame()
    .plotHeatmap(
        object = counts,
        samples = samples,
        genes = genes,
        gene2symbol = gene2symbol(object),
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...)
}



.plotHeatmap.DESeqDataSet <- function(  # nolint
    object,
    normalized = TRUE,
    samples = NULL,
    genes = NULL,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis,
    legendColor = viridis,
    title = NULL,
    ...) {
    .plotHeatmap(
        object = counts(object, normalized = normalized),
        samples = samples,
        genes = genes,
        gene2symbol = gene2symbol,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...)
}



.plotHeatmap.DESeqTransform <- function(  # nolint
    object,
    samples = NULL,
    genes = NULL,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis,
    legendColor = viridis,
    title = NULL,
    ...) {
    .plotHeatmap(
        object = assay(object),
        samples = samples,
        genes = genes,
        gene2symbol = gene2symbol,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...)
}



# Methods ======================================================================
#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("bcbioRNASeq"),
    .plotHeatmap.bcbioRNASeq)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqDataSet"),
    .plotHeatmap.DESeqDataSet)



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqTransform"),
    .plotHeatmap.DESeqTransform)
