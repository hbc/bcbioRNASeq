#' Mapped reads plot
#'
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param bcbio bcbio list object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_mapped_reads(bcbio)
#' }
plot_mapped_reads <- function(bcbio) {
    check_bcbio(bcbio)
    plot <- import_summary(bcbio) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~description,
                          y = ~mapped_reads / 1e6,
                          fill = ~qc_color)
        ) +
        ggplot2::ggtitle("mapped reads") +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(color = "orange",
                            size = 2,
                            yintercept = 10) +
        ggplot2::geom_hline(color = "green",
                            size = 2,
                            yintercept = 20) +
        ggplot2::labs(x = "sample",
                      y = "mapped reads (million)",
                      fill = "") +
        ggplot2::coord_flip()
    print(plot)
}
