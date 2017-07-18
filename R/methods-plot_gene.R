#' Plot an individual gene
#'
#' @rdname plot_gene
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param gene Gene identifier. Can input multiple genes as a character vector.
#' @param format Ensembl identifier format. Defaults to the gene symbol (a.k.a.
#'   `external_gene_name`).
#'
#' @export
setMethod("plot_gene", "bcbioRNADataSet", function(
    object,
    gene,
    format = "symbol") {
    if (!format %in% c("ensgene", "symbol")) {
        stop("Unsupported gene identifier format")
    }
    counts <- tpm(object)
    metadata <- colData(object)
    color <- metadata(object)[["interesting_groups"]][[1L]]

    # Match unique gene identifier with name (gene symbol) using the
    # internally stored Ensembl annotations saved in the run object
    match <- metadata(object)[["annotable"]] %>%
        as.data.frame %>%
        filter(.data[[format]] %in% gene) %>%
        arrange(!!sym("symbol"))

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
            ggtitle(paste(gene_name,
                          paste0("(", gene_id, ")"))) +
            geom_point(size = 4L) +
            theme(
                axis.text.x = element_text(angle = 90L)
            ) +
            labs(x = "sample",
                 y = "transcripts per million (tpm)",
                 fill = "") +
            expand_limits(y = 0L)
        show(plot)
    }) %>% invisible
})

# FIXME Rework to use plot_grid for multiple genes
# FIXME Function needs rework for annotables
# TODO Look at incorporating [plotCounts()] into this function.
