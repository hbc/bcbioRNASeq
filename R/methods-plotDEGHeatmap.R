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
#' @param object Primary object containing [DESeq2::results()] output.
#' @param counts Secondary object containing a normalized counts matrix.
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
#' # Use our stashed gene2symbol for better speed
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
#' @importFrom basejump camel
.plotDEGHeatmap <- function(
    object,
    counts,
    alpha = 0.01,
    lfc = 0,
    title,
    ...) {
    counts <- as.matrix(counts)
    results <- object %>%
        as.data.frame() %>%
        camel(strict = FALSE) %>%
        # Keep genes that pass alpha cutoff
        .[!is.na(.[["padj"]]), , drop = FALSE] %>%
        .[.[["padj"]] < alpha, , drop = FALSE] %>%
        # Keep genes that pass log2 fold change cutoff
        .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
        .[.[["log2FoldChange"]] > lfc |
              .[["log2FoldChange"]] < -lfc, , drop = FALSE]
    if (nrow(results) == 0) {
        warning("No genes passed significance cutoffs", call. = FALSE)
        return(NULL)
    }
    genes <- rownames(results)
    if (length(genes) < 2) {
        warning(paste(
            length(genes), "gene(s) is too few to plot"
        ), call. = FALSE)
        return(NULL)
    } else {
        plotHeatmap(
            counts,
            genes = genes,
            title = title,
            ...)
    }
}



.plotDEGHeatmapDESeqResults <- function(
    object,
    counts,
    lfc = 0,
    title = TRUE,
    ...) {
    results <- as.data.frame(object)
    if (is(counts, "DESeqDataSet") |
        is(counts, "DESeqTransform")) {
        counts <- assay(counts)
    }
    counts <- as.matrix(counts)
    alpha <- metadata(object)[["alpha"]]
    if (isTRUE(title)) {
        title <- .resContrastName(object)
    }
    .plotDEGHeatmap(
        object = results,
        counts = counts,
        alpha = alpha,
        lfc = lfc,
        title = title,
        ...)
}



# Methods ======================================================================
#' @rdname plotDEGHeatmap
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(object = "DESeqResults"),
    .plotDEGHeatmapDESeqResults)
