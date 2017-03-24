#' Plot an individual gene
#'
#' @author Michael Steinbaugh
#'
#' @param gene Gene identifier
#' @param counts Counts matrix
#' @param metadata Metadata data frame
#'
#' @export
plot_gene <- function(gene, counts, metadata) {
    ylab <- deparse(substitute(counts))
    counts <- as.matrix(counts) %>% .[gene, ]
    df <- data.frame(counts = counts,
                     group = metadata$group,
                     name = names(counts))
    plot <- ggplot2::ggplot(
        df,
        ggplot2::aes_(x = ~name,
                      y = ~counts,
                      fill = ~group)
    ) +
        ggplot2::ggtitle(gene) +
        ggplot2::geom_dotplot(binaxis = "y") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90)
        ) +
        ggplot2::xlab("sample") +
        ggplot2::ylab(ylab) +
        ggplot2::expand_limits(y = 0)
    print(plot)
}
