#' Plot an individual gene
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run object
#' @param normalized_counts Normalized counts matrix. Can be obtained from
#'   \code{DESeqDataSet} by running \code{counts(normalized = TRUE)}.
#'   Transcripts per million (TPM) are also acceptable.
#' @param gene Gene identifier
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

    # Convert unique gene identifier to name (gene symbol)
    annotations <- ensembl_annotations(
        run, filters = format, values = gene)
    normalized_counts <- normalized_counts[annotations$ensembl_gene_id, ]
    metadata <- run$metadata

    plot <- data.frame(x = names(normalized_counts),
                       y = normalized_counts,
                       color = metadata[[color]]) %>%
        ggplot(
            aes_(x = ~x,
                 y = ~y,
                 fill = ~color)
        ) +
        ggtitle(gene) +
        geom_dotplot(binaxis = "y") +
        theme(
            axis.text.x = element_text(angle = 90)
        ) +
        labs(x = "sample",
             y = ylab,
             fill = "") +
        expand_limits(y = 0)

    return(plot)
}
