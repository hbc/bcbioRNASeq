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
#' @inherit basejump::plotHeatmap
#'
#' @inheritParams AllGenerics
#' @inheritParams counts
#'
#' @inheritParams gene2symbol
#'
#' @param samples *Optional*. Samples (colnames) to plot.
#' @param genes *Optional*. Gene identifiers (rownames) to plot. These must be
#'   the stable identifiers (e.g. ENSG00000000003) used on Ensembl and not the
#'   gene symbols.
#' @param gene2symbol Apply gene identifier to symbol mappings. If set `TRUE`,
#'   the function will attempt to automatically map gene identifiers to symbols
#'   from Ensembl using [basejump::annotable()]. If set `FALSE`/`NULL`, then
#'   gene2symbol mapping will be disabled. This is useful when working with a
#'   poorly annotated genome. Alternatively, a gene2symbol [data.frame] can be
#'   passed in, and must contain the columns `ensgene` and `symbol`. then the
#'   Ensembl gene identifiers will be labeled in place of gene symbols.
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
    gene2symbol = FALSE,
    annotationCol = NULL,
    scale = "row",
    color = viridis::viridis(256),
    legendColor = viridis::viridis,
    title = NULL,
    ...) {
    # Resize the counts matrix
    if (is.vector(samples)) {
        object <- object[, samples, drop = FALSE]
    }
    if (is.vector(genes)) {
        object <- object[genes, , drop = FALSE]
    }

    # Set the rownames to gene symbols
    if (isTRUE(gene2symbol)) {
        object <- gene2symbol(object)
    } else if (is.data.frame(gene2symbol)) {
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
        title = title)
}


# Methods ======================================================================
#' @rdname plotHeatmap
#' @importFrom viridis viridis
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
        color = viridis::viridis(256),
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
    })



#' @rdname plotHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqDataSet"),
    function(
        object,
        normalized = TRUE,
        samples = NULL,
        genes = NULL,
        gene2symbol = FALSE,
        annotationCol = NULL,
        scale = "row",
        color = viridis::viridis(256),
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
    })



#' @rdname plotHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotHeatmap",
    signature("DESeqTransform"),
    function(
        object,
        genes = NULL,
        gene2symbol = FALSE,
        annotationCol = NULL,
        scale = "row",
        color = viridis::viridis(256),
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
    })
