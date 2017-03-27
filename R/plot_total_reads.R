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
plot_total_reads <- function(summary) {
    summary %>%
        ggplot2::ggplot(aes_(x = ~description,
                             y = ~total_reads / 1e6,
                             fill = ~intgroup)
        ) +
        ggplot2::ggtitle("total reads") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "orange",
                            size = 2,
                            yintercept = 10) +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 20) +
        ggplot2::xlab("sample") +
        ggplot2::ylab("total reads (million)") +
        ggplot2::coord_flip()
}
