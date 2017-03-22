#' Exonic mapping rate plot
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
#' plot_exonic_mapping_rate(summary)
#' }
plot_exonic_mapping_rate <- function(summary) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          # Multiple by 100 here for percentage
                          y = ~exonic_rate * 100,
                          fill = ~group)
        ) +
        ggplot2::ggtitle("exonic mapping rate") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 60) +
        ggplot2::xlab("sample") +
        ggplot2::ylab("exonic mapping rate (%)") +
        ggplot2::ylim(0, 100) +
        ggplot2::coord_flip()
}
