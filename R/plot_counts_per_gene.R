#' @rdname qc_plots
#' @description Counts per gene plot
#'
#' @param run \code{bcbio-nextgen} run object
#' @param normalized_counts Normalized counts matrix
#'
#' @return Boxplot
#' @export
plot_counts_per_gene <- function(run, normalized_counts) {
    check_run(run)

    color <- run$intgroup[1]
    name <- deparse(substitute(normalized_counts))

    melted <- melt_log10(run, normalized_counts)

    plot <- data.frame(x = melted$description,
                       y = melted$counts,
                       color = melted[[color]]) %>%
        ggplot(
            aes_(x = ~x,
                 y = ~y,
                 color = ~color)
        ) +
        ggtitle(paste("counts per gene:", name)) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "sample",
             y = expression(log[10]~counts~per~gene),
             color = "") +
        coord_flip()

    return(plot)
}
