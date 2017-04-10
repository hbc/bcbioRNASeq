#' @rdname qc_plots
#' @description Counts per gene plot
#' @return Boxplot
#' @export
plot_counts_per_gene <- function(bcbio, normalized_counts) {
    check_bcbio(bcbio)

    color <- bcbio$intgroup[1]
    name <- deparse(substitute(normalized_counts))

    melted <- melt_log10(bcbio, normalized_counts)

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
