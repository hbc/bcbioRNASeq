#' Plot an individual gene
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run object
#' @param normalized_counts Normalized counts matrix. Can be obtained from
#'   \code{DESeqDataSet} by running \code{counts(normalized = TRUE)}.
#'   Transcripts per million (TPM) are also acceptable.
#' @param gene Gene identifier
#'
#' @export
plot_gene <- function(
    run,
    normalized_counts,
    gene) {
    check_run(run)

    color <- run$intgroup[1]
    ylab <- deparse(substitute(normalized_counts))

    normalized_counts <- as.matrix(normalized_counts) %>% .[gene, ]
    metadata <- import_metadata(run)

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
