#' Count density plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param bcbio bcbio list object
#' @param counts Counts matrix
#' @param print Print plot
#'
#' @return Density plot
#' @export
#'
#' @examples
#' \dontrun{
#' plot_count_density(bcbio, counts)
#' }
plot_count_density <- function(
    bcbio,
    counts,
    print = TRUE) {
    name <- deparse(substitute(counts))
    melted <- melt_log10(bcbio, counts)
    plot <- ggplot2::ggplot(
        melted,
        ggplot2::aes_(x = ~counts,
                      # Use `~description` instead?
                      group = ~samplename)
    ) +
        ggplot2::ggtitle(paste("count density:", name)) +
        ggplot2::geom_density() +
        ggplot2::labs(x = expression(log[10]~counts~per~gene),
                      y = "density") +
    if (isTRUE(print)) {
        print(plot)
    } else {
        return(plot)
    }
}
