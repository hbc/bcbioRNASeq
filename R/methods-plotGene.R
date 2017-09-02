#' Plot Individual Genes
#'
#' @rdname plotGene
#' @name plotGene
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
#' plotGene(bcb, gene = c("Sulf1", "Phf3"))
NULL



# Methods ====
#' @rdname plotGene
#' @export
setMethod("plotGene", "bcbioRNADataSet", function(
    object,
    normalized = "tpm",
    gene,
    format = "symbol") {
    if (!format %in% c("ensgene", "symbol")) {
        stop("Unsupported gene identifier format")
    }
    if (is.logical(normalized)) {
        stop("Explicit request of normalized counts format is required")
    }
    counts <- counts(object, normalized = normalized)
    metadata <- colData(object)
    color <- metadata(object)[["interestingGroups"]][[1L]]

    # Match unique gene identifier with name (gene symbol) using the
    # internally stored Ensembl annotations saved in the run object
    match <- metadata(object)[["annotable"]] %>%
        .[.[[format]] %in% gene, ] %>%
        .[, c("symbol", "ensgene")]

    # Seq along Ensembl data frame here instead of the gene input vector,
    # which will then output only genes that match Ensembl
    lapply(1L:nrow(match), function(a) {
        geneName <- match[["symbol"]][[a]]
        geneID <- match[["ensgene"]][[a]]
        plot <- data.frame(
            x = colnames(counts),
            y = counts[geneID, ],
            color = metadata[[color]]) %>%
            ggplot(
                aes_(x = ~x,
                     y = ~y,
                     color = ~color)
            ) +
            geom_point(size = 4L) +
            theme(
                axis.text.x = element_text(angle = 90L)
            ) +
            labs(title = geneName,
                 x = "sample",
                 y = paste0("counts (", normalized, ")"),
                 color = color) +
            expand_limits(y = 0L) +
            scale_color_viridis(discrete = TRUE)
        show(plot)
    }) %>%
        invisible
})
