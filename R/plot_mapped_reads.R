#' Mapped reads plot
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
#' plot_mapped_reads(summary)
#' }
plot_mapped_reads <- function(summary) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = ~mapped_reads / 1e6,
                          fill = ~intgroup)
        ) +
        ggplot2::ggtitle("mapped reads") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "orange",
                            size = 2,
                            yintercept = 10) +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 20) +
        ggplot2::xlab("sample") +
        ggplot2::ylab("mapped reads (million)") +
        ggplot2::coord_flip()
}
