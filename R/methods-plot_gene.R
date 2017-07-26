#' Plot Individual Genes
#'
#' @rdname plot_gene
#' @author Michael Steinbaugh
#'
#' @param gene Gene identifier. Can input multiple genes as a character vector.
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @param format Ensembl identifier format. Supports `symbol` (**default**) or
#'   `ensgene`.
#'
#' @return [ggplot].
#' @export
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' data(bcb)
#' plot_gene(bcb, gene = c("Sulf1", "Phf3"))
setMethod("plot_gene", "bcbioRNADataSet", function(
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
    color <- metadata(object)[["interesting_groups"]][[1L]]

    # Match unique gene identifier with name (gene symbol) using the
    # internally stored Ensembl annotations saved in the run object
    match <- metadata(object)[["annotable"]] %>%
        filter(.data[[format]] %in% !!gene) %>%
        tidy_select(c("symbol", "ensgene"))

    # Seq along Ensembl data frame here instead of the gene input vector,
    # which will then output only genes that match Ensembl
    lapply(1L:nrow(match), function(a) {
        gene_name <- match[["symbol"]][[a]]
        gene_id <- match[["ensgene"]][[a]]
        plot <- data.frame(
            x = colnames(counts),
            y = counts[gene_id, ],
            color = metadata[[color]]) %>%
            ggplot(
                aes_(x = ~x,
                     y = ~y,
                     color = ~color,
                     shape = ~color)
            ) +
            geom_point(size = 4L) +
            theme(
                axis.text.x = element_text(angle = 90L)
            ) +
            labs(title = gene_name,
                 x = "sample",
                 y = str_c("counts (", normalized, ")"),
                 color = color,
                 shape = color) +
            expand_limits(y = 0L)
        show(plot)
    }) %>% invisible
})
