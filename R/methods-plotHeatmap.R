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
#' @importFrom bcbioBase plotHeatmap
#'
#' @inheritParams general
#' @inheritParams counts
#' @inheritParams gene2symbol
#' @inheritParams plotGene
#'
#' @param samples *Optional*. Samples (colnames) to plot.
#' @param annotationCol *Optional*. [data.frame] that defines annotation
#'   mappings for the columns.
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the row direction or the column direction, or none. Corresponding
#'   values are "row", "column" and "none".
#' @param legendColor Colors to use for legend labels. Defaults to the
#'   [viridis::viridis()].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # Use Ensembl identifiers to define genes
#' ensgene <- rownames(bcb)[1:20]
#' plotHeatmap(bcb, genes = ensgene)
#'
#' # Use inferno color palette
#' plotHeatmap(
#'     bcb,
#'     genes = ensgene,
#'     color = viridis::inferno(256),
#'     legendColor = viridis::inferno)
#'
#' # Transcriptome heatmap with default pheatmap colors
#' plotHeatmap(
#'     bcb,
#'     color = NULL,
#'     legendColor = NULL)
#'
#' # DESeqDataSet
#' dds <- bcbio(bcb, "DESeqDataSet")
#' plotHeatmap(dds)
#'
#' # DESeqTransform
#' rld <- assays(bcb)[["rlog"]]
#' plotHeatmap(rld)
NULL



# Constructors =================================================================
.plotHeatmap <- function(
    object,
    samples = NULL,
    genes = NULL,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis::viridis,
    legendColor = viridis::viridis,
    title = NULL,
    ...) {
    assert_is_matrix(object)
    # Passthrough: color, legendColor
    .assert_gene2symbol(object, genes, gene2symbol)
    # TODO Migrate to `assert_formal_color_function`
    assert_is_any_of(color, c("function", "NULL"))
    assert_is_any_of(legendColor, c("function", "NULL"))

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

    # TODO Safe to remove once basejump is updated
    if (is.function(color)) {
        color <- color(256L)
    }

    plotHeatmap(
        object = object,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title)
}



#' @importFrom viridis viridis
.plotHeatmap.bcbioRNASeq <- function(  # nolint
    object,
    normalized = "rlog",
    samples = NULL,
    genes = NULL,
    scale = "row",
    color = viridis::viridis,
    legendColor = viridis::viridis,
    title = NULL,
    ...) {
    counts <- counts(object, normalized = normalized)
    gene2symbol <- gene2symbol(object)
    annotationCol <- colData(object) %>%
        .[colnames(counts), interestingGroups(object), drop = FALSE]
    .plotHeatmap(
        object = counts,
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



#' @importFrom viridis viridis
.plotHeatmap.DESeqDataSet <- function(
    object,
    normalized = TRUE,
    samples = NULL,
    genes = NULL,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis::viridis,
    legendColor = viridis::viridis,
    title = NULL,
    ...) {
    counts <- counts(object, normalized = normalized)
    .plotHeatmap(
        object = counts,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...)
}



#' @importFrom viridis viridis
.plotHeatmap.DESeqTransform <- function(  # nolint
    object,
    genes = NULL,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis::viridis(256L),
    legendColor = viridis::viridis,
    title = NULL,
    ...) {
    counts <- assay(object)
    .plotHeatmap(
        object = counts,
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
