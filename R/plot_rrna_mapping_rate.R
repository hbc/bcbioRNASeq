#' rRNA contamination mapping rate
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
#' plot_rrna_mapping_rate(summary)
#' }
plot_rrna_mapping_rate <- function(summary) {
    summary %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = ~rrna_rate * 100,
                          fill = ~group)) +
        ggplot2::ggtitle("rRNA mapping rate") +
        ggplot2::geom_bar(stat = "identity") +
        #` ggplot2::geom_hline(linetype = 2, yintercept = 5) +
        ggplot2::geom_hline(color = "red",
                            size = 2,
                            yintercept = 10) +
        ggplot2::xlab("sample") +
        ggplot2::ylab("rRNA mapping rate (%)") +
        #` ggplot2::ylim(0, 100) +
        ggplot2::coord_flip()
}
