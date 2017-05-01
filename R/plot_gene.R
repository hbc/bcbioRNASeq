#' Plot an individual gene
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run object
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
    ensembl <- run$ensembl %>% .[.[[format]] %in% gene, ]
    normalized_counts <- normalized_counts[ensembl$ensembl_gene_id, ]

    # Seq along Ensembl data frame here instead of the gene input vector,
    # which will then output only genes that match Ensembl
    lapply(seq_along(nrow(ensembl)), function(a) {
        plot <- data.frame(x = colnames(normalized_counts),
                           y = normalized_counts[a, ],
                           color = metadata[[color]]) %>%
            ggplot(
                aes_(x = ~x,
                     y = ~y,
                     fill = ~color)
            ) +
            ggtitle(gene[a]) +
            geom_dotplot(binaxis = "y") +
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
