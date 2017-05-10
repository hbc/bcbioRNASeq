#' Plot an individual gene
#'
#' @author Michael Steinbaugh
#'
#' @param run bcbio-nextgen run object
#' @param normalized_counts Normalized counts matrix. Can be obtained from
#'   \code{DESeqDataSet} by running \code{counts(normalized = TRUE)}.
#'   Transcripts per million (TPM) are also acceptable.
#' @param gene Gene identifier. Can input multiple genes as a character vector.
#' @param format Ensembl identifier format. Defaults to the gene name (symbol).
#'
#' @export
plot_gene <- function(
    run,
    normalized_counts,
    gene,
    format = "external_gene_name") {
    check_run(run)

    color <- run$intgroup[1]
    ylab <- deparse(substitute(normalized_counts))
    metadata <- run$metadata

    # Match unique gene identifier with name (gene symbol) using the internally
    # stored Ensembl annotations saved in the run object
    match <- run$ensembl %>%
        .[.[[format]] %in% gene, ] %>%
        .[order(.$external_gene_name), ]

    # Seq along Ensembl data frame here instead of the gene input vector,
    # which will then output only genes that match Ensembl
    lapply(1:nrow(match), function(a) {
        gene_name <- match$external_gene_name[[a]]
        gene_id <- match$ensembl_gene_id[[a]]
        plot <- data.frame(
            x = colnames(normalized_counts),
            y = normalized_counts[gene_id, ],
            color = metadata[[color]]) %>%
            ggplot(
                aes_(x = ~x,
                     y = ~y,
                     color = ~color,
                     shape = ~color)
            ) +
            ggtitle(paste(gene_name, gene_id, sep = " : ")) +
            geom_point(size = 4) +
            theme(
                axis.text.x = element_text(angle = 90)
            ) +
            labs(x = "sample",
                 y = ylab,
                 fill = "") +
            expand_limits(y = 0)
        show(plot)
    }) %>% invisible
}
