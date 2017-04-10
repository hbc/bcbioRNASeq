#' Plot an individual gene
#'
#' @author Michael Steinbaugh
#'
#' @param bcbio bcbio list object
#' @param counts Counts matrix
#' @param gene Gene identifier
#'
#' @export
plot_gene <- function(
    bcbio,
    counts,
    gene) {
    check_bcbio(bcbio)
    color <- bcbio$intgroup[1]
    ylab <- deparse(substitute(counts))
    counts <- as.matrix(counts) %>% .[gene, ]

    # Import metadata with automatic lane split detection. This will match the
    # colnames format of the input counts.
    if (any(grepl("_L\\d+$", colnames(counts)))) {
        lanes <- "split"
    } else {
        lanes <- NULL
    }
    metadata <- import_metadata(bcbio)

    plot <- data.frame(x = names(counts),
                       y = counts,
                       color = metadata[[color]]) %>%
        ggplot(
            aes_(x = ~name,
                 y = ~counts,
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
