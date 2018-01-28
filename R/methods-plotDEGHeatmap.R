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
#' @param counts Secondary object containing a normalized count matrix.
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
#'
#' # DESeqResults, DESeqTransform
#' plotDEGHeatmap(
#'     res,
#'     counts = rld,
#'     gene2symbol = gene2symbol)
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
    results,
    counts,
    alpha = 0.01,
    lfc = 0,
    title,
    ...) {
    results <- results %>%
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



# Methods ======================================================================
#' @rdname plotDEGHeatmap
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        object = "DESeqResults",
        counts = "DESeqTransform"),
    function(
        object,
        counts,
        lfc = 0,
        title = TRUE,
        ...) {
        results <- as.data.frame(object)
        counts <- assay(counts)
        alpha <- metadata(object)[["alpha"]]
        if (isTRUE(title)) {
            title <- .resContrastName(object)
        }
        .plotDEGHeatmap(
            results = results,
            counts = counts,
            alpha = alpha,
            lfc = lfc,
            title = title,
            ...)
    })



#' @rdname plotDEGHeatmap
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(
        object = "DESeqResults",
        counts = "DESeqDataSet"),
    function(
        object,
        counts,
        lfc = 0,
        title = TRUE,
        ...) {
        warning("DESeqTransform for counts is recommended", call. = FALSE)
        results <- as.data.frame(object)
        counts <- counts(counts, normalized = TRUE)
        alpha <- metadata(object)[["alpha"]]
        if (isTRUE(title)) {
            title <- .resContrastName(object)
        }
        .plotDEGHeatmap(
            results = results,
            counts = counts,
            alpha = alpha,
            lfc = lfc,
            title = title,
            ...)
    })
