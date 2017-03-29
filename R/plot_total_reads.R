#' Total reads plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param summary \code{bcbio-rnaseq} summary report
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_total_reads(summary)
#' }
plot_total_reads <- function(bcbio, summary) {
    color <- bcbio$intgroup[1]
    data.frame(x = summary$description,
               y = summary$total_reads / 1e6,
               fill = summary[[color]]) %>%
        ggplot2::ggplot(aes_(x = ~x,
                             y = ~y,
                             fill = ~fill)
        ) +
        ggplot2::ggtitle("total reads") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "orange",
                            size = 2,
                            yintercept = 10) +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 20) +
        ggplot2::labs(x = "sample",
                      y = "total reads (million)",
                      fill = "") +
        ggplot2::coord_flip()
}
