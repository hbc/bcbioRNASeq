#' Plot Heatmap
#'
#' @details
#' When column clustering is disabled, the columns are sorted by the interesting
#' groups (`interestingGroups`) specified in the [bcbioRNASeq] and then the
#' sample names.
#'
#' @name plotHeatmap
#' @family Gene Expression Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotHeatmap
#'
#' @inherit basejump::plotHeatmap
#'
#' @inheritParams general
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/rld_small.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq ====
#' plotHeatmap(bcb_small, genes = head(rownames(bcb_small), 20L))
#'
#' # Full transcriptome heatmap with default pheatmap colors
#' plotHeatmap(bcb_small, color = inferno, legendColor = inferno)
#'
#' # DESeqDataSet ====
#' plotHeatmap(
#'     object = dds_small,
#'     genes = head(rownames(dds_small), 20L),
#'     gene2symbol = gene2symbol(bcb_small)
#' )
#'
#' # DESeqTransform ====
#' plotHeatmap(
#'     object = rld_small,
#'     genes = head(rownames(rld_small), 20L),
#'     gene2symbol = gene2symbol(bcb_small)
#' )
NULL



# Constructors =================================================================
.plotHeatmap.matrix <- function(  # nolint
    object,
    samples = NULL,
    genes = NULL,
    gene2symbol = NULL,
    annotationCol = NULL,
    scale = "row",
    color = viridis,
    legendColor = viridis,
    title = NULL,
    ...
) {
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
            gene2symbol = gene2symbol
        )
    }

    plotHeatmap(
        object = object,
        annotationCol = annotationCol,
        scale = scale,
        color = color,
        legendColor = legendColor,
        title = title,
        ...
    )
}



# Methods ======================================================================
#' @rdname plotHeatmap
#' @export
setMethod(
    "plotHeatmap",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("rlog", "vst", "tmm", "tpm"),
        samples = NULL,
        genes = NULL,
        scale = "row",
        color = viridis,
        legendColor = viridis,
        title = NULL,
        ...
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)
        annotationCol <- colData(object) %>%
            .[colnames(counts), interestingGroups(object), drop = FALSE] %>%
            as.data.frame()
        .plotHeatmap.matrix(
            object = counts,
            samples = samples,
            genes = genes,
            gene2symbol = gene2symbol(object),
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            ...
        )
    }
)



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
        ...
    ) {
        validObject(object)
        assert_is_a_bool(normalized)
        .plotHeatmap.matrix(
            object = counts(object, normalized = normalized),
            samples = samples,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            ...
        )
    }
)



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
        ...
    ) {
        validObject(object)
        .plotHeatmap.matrix(
            object = assay(object),
            samples = samples,
            genes = genes,
            gene2symbol = gene2symbol,
            annotationCol = annotationCol,
            scale = scale,
            color = color,
            legendColor = legendColor,
            title = title,
            ...
        )
    }
)
