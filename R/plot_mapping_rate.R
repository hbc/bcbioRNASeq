#' Mapping rate plot
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
#' plot_mapping_rate(summary)
#' }
plot_mapping_rate <- function(summary) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = ~mapped_reads / total_reads * 100,
                          fill = ~intgroup)
        ) +
        ggplot2::ggtitle("mapping rate") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "orange",
                            size = 2,
                            yintercept = 70) +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 90) +
        ggplot2::xlab("sample") +
        ggplot2::ylab("mapping rate (%)") +
        ggplot2::ylim(0, 100) +
        ggplot2::coord_flip()
}
