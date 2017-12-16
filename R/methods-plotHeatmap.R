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
#' @inheritParams AllGenerics
#' @inheritParams basejump::plotHeatmap
#' @inheritParams counts
#'
#' @return [pheatmap::pheatmap()] ggplot return.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # Use Ensembl identifiers to define genes
#' ensgene <- counts(bcb)[1:20, ] %>% rownames()
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
        genes = NULL,
        scale = "row",
        color = viridis::viridis(256),
        legendColor = viridis::viridis,
        title = NULL,
        quiet = FALSE,
        ...) {
        counts <- counts(object, normalized = normalized)
        gene2symbol <- gene2symbol(object)
        interestingGroups <- interestingGroups(object)
        annotationCol <- colData(object) %>%
            .[, interestingGroups, drop = FALSE]
        plotHeatmap(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            quiet = quiet,
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
        genes = NULL,
        gene2symbol = TRUE,
        annotationCol = NULL,
        scale = "row",
        color = viridis::viridis(256),
        legendColor = viridis::viridis,
        title = NULL,
        quiet = FALSE,
        ...) {
        counts <- counts(object, normalized = normalized)
        plotHeatmap(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            quiet = quiet,
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
        gene2symbol = TRUE,
        annotationCol = NULL,
        scale = "row",
        color = viridis::viridis(256),
        legendColor = viridis::viridis,
        title = NULL,
        quiet = FALSE,
        ...) {
        counts <- assay(object)
        plotHeatmap(
            object = counts,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            quiet = quiet,
            ...)
    })
