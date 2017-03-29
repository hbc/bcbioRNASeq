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
    color <- bcbio$intgroup[1]
    ylab <- deparse(substitute(counts))
    counts <- as.matrix(counts) %>% .[gene, ]

    # Import metadata with automatic lane split detection. This will match the
    # colnames format of the input counts.
    if (any(grepl("_L\\d+$", colnames(counts)))) {
        lane_split <- TRUE
    } else {
        lane_split <- FALSE
    }
    metadata <- import_metadata(bcbio, lane_split = lane_split)

    plot <- data.frame(x = names(counts),
                       y = counts,
                       color = metadata[[color]]) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~name,
                          y = ~counts,
                          fill = ~color)
        ) +
        ggplot2::ggtitle(gene) +
        ggplot2::geom_dotplot(binaxis = "y") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90)
        ) +
        ggplot2::labs(x = "sample",
                      y = ylab,
                      fill = "") +
        ggplot2::expand_limits(y = 0)
    print(plot)
}
