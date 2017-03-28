#' Counts per gene plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param bcbio bcbio list object
#' @param counts Counts matrix
#' @param flip Flip x and y axes
#' @param print Print plot
#'
#' @return Boxplot
#' @export
#'
#' @examples
#' \dontrun{
#' plot_counts_per_gene(bcbio, counts)
#' }
plot_counts_per_gene <- function(
    bcbio,
    counts,
    flip = TRUE,
    print = TRUE) {
    name <- deparse(substitute(counts))
    melted <- melt_log10(bcbio, counts)
    plot <- data.frame(x = melted$samplename,
                       y = melted$counts,
                       color = melted[[bcbio$intgroup[1]]]) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~x,
                          y = ~y,
                          color = ~color)
        ) +
        ggplot2::ggtitle(paste("counts per gene:", name)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::xlab("sample") +
        ggplot2::ylab(expression(log[10]~counts~per~gene))
    if (isTRUE(flip)) {
        plot <- plot +
            ggplot2::coord_flip()
    }
    if (isTRUE(print)) {
        print(plot)
    } else {
        return(plot)
    }
}
