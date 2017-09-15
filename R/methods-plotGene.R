#' Plot Individual Genes
#'
#' @rdname plotGene
#' @name plotGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @param gene Gene identifier. Can input multiple genes as a character vector.
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @param format Ensembl identifier format. Supports `symbol` (**default**) or
#'   `ensgene`.
#'
#' @return [ggplot].
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' data(bcb)
#' genes <- c("Sulf1", "Phf3")
#'
#' # bcbioRNADataSet
#' plotGene(bcb, gene = genes)
#' plotGene(bcb, gene = genes, interestingGroup = "group")
NULL



# Constructors ====
.plotGene <- function(
    counts,
    gene,
    metadata,
    interestingGroup = "sampleName") {
    metadata <- as.data.frame(metadata)
    sapply(seq_along(gene), function(a) {
        ensgene <- gene[[a]]
        symbol <- names(gene)[[a]]
        df <- data.frame(
            x = colnames(counts),
            y = counts[ensgene, ],
            color = metadata[[interestingGroup]])
        p <- ggplot(df,
                aes_(x = ~x,
                     y = ~y,
                     color = ~color)) +
            geom_point(size = 4L) +
            theme(
                axis.text.x = element_text(angle = 90L)) +
            labs(title = symbol,
                 x = "sample",
                 y = paste0("counts (", normalized, ")"),
                 color = interestingGroup) +
            expand_limits(y = 0L) +
            scale_color_viridis(discrete = TRUE)
        show(p)
    }) %>%
        invisible
}



# Methods ====
#' @rdname plotGene
#' @export
setMethod("plotGene", "bcbioRNADataSet", function(
    object,
    interestingGroup,
    normalized = "tpm",
    gene,
    format = "symbol") {
    if (!format %in% c("ensgene", "symbol")) {
        stop("Unsupported gene identifier format")
    }
    if (is.logical(normalized)) {
        warning("Explicit format of normalized counts format is recommended")
    }
    if (missing(interestingGroup)) {
        interestingGroup <- .interestingGroup(object)
    }

    counts <- counts(object, normalized = normalized)
    metadata <- colData(object)

    # Match unique gene identifier with name (gene symbol) using the
    # internally stored Ensembl annotations saved in the run object
    gene2symbol <- metadata(object) %>%
        .[["annotable"]] %>%
        .[.[[format]] %in% gene, , drop = FALSE] %>%
        .[, c("symbol", "ensgene")]
    gene <- gene2symbol[["ensgene"]]
    names(gene) <- gene2symbol[["symbol"]]

    .plotGene(
        counts = counts,
        gene = gene,
        metadata = metadata,
        interestingGroup = interestingGroup)
})
