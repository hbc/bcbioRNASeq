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
#'     object = dds,
#'     genes = head(rownames(dds), 20L),
#'     gene2symbol = gene2symbol)
#'
#' # DESeqTransform ====
#' rld <- assays(bcb)[["rlog"]]
#' plotHeatmap(
#'     object = rld,
#'     genes = head(rownames(rld), 20L),
#'     gene2symbol = gene2symbol)
NULL



# Constructors =================================================================
#' @importFrom basejump convertGenesToSymbols plotHeatmap
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
    assertIsCharacterOrNULL(samples)
    assertIsCharacterOrNULL(genes)
    assertFormalGene2symbol(object, genes, gene2symbol)
    assertFormalAnnotationCol(object, annotationCol)
    assertIsHexColorFunctionOrNULL(color)
    assertIsHexColorFunctionOrNULL(legendColor)
    assertIsAStringOrNULL(title)

    # Resize the counts matrix
    if (is.vector(samples)) {
        object <- object[, samples, drop = FALSE]
    }
    if (is.vector(genes)) {
        object <- object[genes, , drop = FALSE]
    }

    # Set the rownames to gene symbols
    if (is.data.frame(gene2symbol)) {
        rownames(object) <- convertGenesToSymbols(
            rownames(object),
            gene2symbol = gene2symbol)
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



# Methods ======================================================================
#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("bcbioRNASeq"),
    function(
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
    })



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqDataSet"),
    function(
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
    })



#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqTransform"),
    function(
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
    })
