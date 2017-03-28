#' Count density plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param counts Counts matrix
#' @param metadata Metadata data frame
#' @param print Print plot
#'
#' @return Density plot
#' @export
#'
#' @examples
#' \dontrun{
#' plot_count_density(counts, metadata)
#' }
plot_count_density <- function(
    counts,
    metadata,
    print = TRUE) {
    name <- deparse(substitute(counts))
    plot <- counts %>%
        melt_log10(metadata = metadata) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~counts,
                          # Use `~description` instead?
                          group = ~samplename)
        ) +
        ggplot2::ggtitle(paste("count density:", name)) +
        ggplot2::geom_density() +
        ggplot2::xlab(expression(log[10]~counts~per~gene)) +
        ggplot2::ylab("density")
    if (isTRUE(print)) {
        print(plot)
    } else {
        return(plot)
    }
}
